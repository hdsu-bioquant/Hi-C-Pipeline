# =============================================================================
# align.smk — bwa-mem alignment streamed straight into pairtools parse + sort
#
# The 4D Nucleome standard for Hi-C:
#   bwa mem -SP5M  (align the two read ends independently; 5'-primary for the
#                   chimeric reads created by the ligation junction)
#   | pairtools parse  (classify each read pair -> a .pairs record + pair_type)
#   | pairtools sort   (block-sort by (chr1,chr2,pos1,pos2) for downstream binning)
#
# Streaming avoids writing a multi-hundred-GB intermediate BAM.
# Alignment is per (sample, lane-unit); lanes are merged in pairs.smk before dedup.
# =============================================================================

rule align_parse_sort:
    input:
        unpack(fq_for_align),
        idx = multiext(FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        cs  = CHROMSIZES
    output:
        pairs = temp(str(RESULTS / "pairs/{sample}/{unit}.sorted.pairs.gz"))
    params:
        fa            = FASTA,
        bwa_flags     = config["bwa"]["flags"],
        walks         = config["pairtools"]["walks_policy"],
        max_gap       = config["pairtools"]["max_inter_align_gap"],
        min_mapq      = config["pairtools"]["min_mapq"],
        sort_threads  = config["pairtools"]["sort_threads"],
        sort_mem      = config["pairtools"]["sort_mem"],
        tmpdir        = lambda wc: str(RESULTS / "pairs" / f"tmp_{wc.sample}_{wc.unit}")
    log:
        str(LOGS / "align/{sample}.{unit}.bwa_pairtools.log")
    container:
        IMG_4DN
    threads: config["bwa"]["threads"]
    resources:
        mem_mb = 64000, runtime = 2880
    shell:
        r"""
        mkdir -p $(dirname {output.pairs}) {params.tmpdir}
        ( bwa mem {params.bwa_flags} -t {threads} {params.fa} {input.r1} {input.r2} \
          | pairtools parse \
              --chroms-path {input.cs} \
              --assembly {GENOME} \
              --min-mapq {params.min_mapq} \
              --walks-policy {params.walks} \
              --max-inter-align-gap {params.max_gap} \
              --add-columns mapq \
              --drop-sam \
          | pairtools sort \
              --nproc {params.sort_threads} \
              --memory {params.sort_mem} \
              --tmpdir {params.tmpdir} \
              -o {output.pairs} \
        ) > {log} 2>&1
        rm -rf {params.tmpdir}
        """
