# =============================================================================
# pairs.smk — lane merge, PCR/optical dedup, valid-pair selection, pairix index
#
#   per-lane sorted.pairs -> merge_lanes (pairtools merge, keeps sort order)
#                         -> dedup (remove duplicates, emit stats)
#                         -> select valid pair types (UU/UR/RU)
#                         -> bgzip + pairix index  (ready for cooler cload pairix)
#
# pairtools dedup stats are MultiQC-readable and give you the key library
# metrics: % duplicates, cis/trans ratio, valid-pair yield.
# =============================================================================

rule merge_lanes:
    """
    Merge the per-lane sorted .pairs of one library into a single sorted stream.
    `pairtools merge` is a sorted merge, so the result stays correctly ordered for
    dedup. For a single-lane library this is effectively a pass-through.
    """
    input:
        lanes = units_sorted_pairs
    output:
        pairs = temp(str(RESULTS / "pairs/{sample}.sorted.pairs.gz"))
    log:
        str(LOGS / "pairs/{sample}.merge_lanes.log")
    container:
        IMG_4DN
    threads: 4
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        pairtools merge --nproc {threads} -o {output.pairs} {input.lanes} > {log} 2>&1
        """


rule pairs_dedup:
    input:
        pairs = str(RESULTS / "pairs/{sample}.sorted.pairs.gz")
    output:
        pairs = temp(str(RESULTS / "pairs/{sample}.dedup.pairs.gz")),
        stats = str(RESULTS / "pairs/{sample}.dedup.stats.txt")
    log:
        str(LOGS / "pairs/{sample}.dedup.log")
    container:
        IMG_4DN
    threads: 4
    resources:
        mem_mb = 16000, runtime = 480
    shell:
        r"""
        pairtools dedup \
            --mark-dups \
            --output-stats {output.stats} \
            --output {output.pairs} \
            {input.pairs} > {log} 2>&1
        """


rule pairs_select:
    input:
        pairs = str(RESULTS / "pairs/{sample}.dedup.pairs.gz")
    output:
        pairs = str(RESULTS / "pairs/{sample}.valid.pairs.gz"),
        index = str(RESULTS / "pairs/{sample}.valid.pairs.gz.px2")
    params:
        expr = config["pairtools"]["select_expr"]
    log:
        str(LOGS / "pairs/{sample}.select.log")
    container:
        IMG_4DN
    threads: 4
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        # NOTE: the duplexa/4dn-hic:v43 image ships an OLD bgzip (htslib < 1.4)
        # that has NO threads flag — `bgzip -@ N` fails with
        # "bgzip: invalid option -- '@'". Use plain single-threaded bgzip; it
        # still produces BGZF that pairix can index. (Don't re-add `-@` unless
        # this rule moves to a newer htslib image.)
        ( pairtools select '{params.expr}' {input.pairs} \
            | bgzip > {output.pairs}
          pairix -f {output.pairs}
        ) > {log} 2>&1
        """
