# =============================================================================
# qc.smk — read QC (FastQC), optional trimming (fastp), aggregation (MultiQC)
#
# FastQC/fastp run per (sample, lane-unit) so multi-lane libraries are QC'd
# correctly. MultiQC then aggregates everything (FastQC + pairtools dedup stats).
#
# NOTE on trimming: Hi-C generally needs NO read trimming. The ligation junction
# produces chimeric reads that are handled at the alignment/pairing stage by
# `bwa mem -SP5M` + `pairtools --walks-policy 5unique`, not by trimming.
# `trim` is therefore OFF by default; enable only for adapter contamination.
# =============================================================================

rule fastqc:
    input:
        fq = lambda wc: unit_fastq(wc.sample, wc.unit, "r1" if wc.r == "R1" else "r2")
    output:
        zip  = str(RESULTS / "qc/fastqc/{sample}.{unit}_{r}_fastqc.zip"),
        html = str(RESULTS / "qc/fastqc/{sample}.{unit}_{r}_fastqc.html")
    log:
        str(LOGS / "qc/fastqc_{sample}.{unit}_{r}.log")
    container:
        IMG_FASTQC
    threads: 2
    resources:
        mem_mb = 4000, runtime = 120
    shell:
        r"""
        mkdir -p $(dirname {output.zip})
        # FastQC names outputs after the input basename; rename to sample/unit/read.
        tmp=$(mktemp -d)
        fastqc -t {threads} -o "$tmp" {input.fq} > {log} 2>&1
        mv "$tmp"/*_fastqc.zip  {output.zip}
        mv "$tmp"/*_fastqc.html {output.html}
        rm -rf "$tmp"
        """


rule fastp_trim:
    """Optional adapter-only trimming (only used when runrule.trim = true)."""
    input:
        r1 = lambda wc: unit_fastq(wc.sample, wc.unit, "r1"),
        r2 = lambda wc: unit_fastq(wc.sample, wc.unit, "r2")
    output:
        r1   = str(RESULTS / "trim/{sample}.{unit}_R1.trimmed.fastq.gz"),
        r2   = str(RESULTS / "trim/{sample}.{unit}_R2.trimmed.fastq.gz"),
        json = str(RESULTS / "trim/{sample}.{unit}.fastp.json"),
        html = str(RESULTS / "trim/{sample}.{unit}.fastp.html")
    log:
        str(LOGS / "qc/fastp_{sample}.{unit}.log")
    container:
        IMG_FASTP
    threads: 8
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        mkdir -p $(dirname {output.r1})
        # Adapter trimming only; do NOT quality-crop Hi-C reads aggressively.
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --disable_quality_filtering --disable_length_filtering \
              --thread {threads} -j {output.json} -h {output.html} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc = [str(RESULTS/f"qc/fastqc/{s}.{u}_{r}_fastqc.zip")
                  for (s, u) in UNITS for r in ["R1", "R2"]] if run("fastqc") else [],
        stats  = expand(str(RESULTS/"pairs/{sample}.dedup.stats.txt"),
                        sample=SAMPLES) if run("align_pairs") else []
    output:
        html = str(RESULTS / "qc/multiqc/multiqc_report.html")
    params:
        fastqc_dir = str(RESULTS / "qc/fastqc"),
        pairs_dir  = str(RESULTS / "pairs")
    log:
        str(LOGS / "qc/multiqc.log")
    container:
        IMG_MULTIQC
    threads: 1
    resources:
        mem_mb = 4000, runtime = 60
    shell:
        r"""
        mkdir -p $(dirname {output.html})
        # MultiQC understands FastQC and pairtools stats natively.
        multiqc -f -o $(dirname {output.html}) -n multiqc_report.html \
            {params.fastqc_dir} {params.pairs_dir} > {log} 2>&1
        """
