# =============================================================================
# matrix.smk — build contact matrices in BOTH community-standard formats
#
#   .mcool  (cooler)  -> for cooltools (insulation TADs, eigs compartments, saddle)
#                        ICE-balanced multi-resolution pyramid.
#   .hic    (Juicer)  -> for Juicer tools (Arrowhead, HiCCUPS, APA, HiCCUPSDiff)
#                        KR + SCALE normalisations computed by juicer_tools pre.
#
# Both are derived from the SAME validated .pairs file, so they describe identical
# contacts -- they differ only in container format and balancing convention.
# =============================================================================

rule cooler_cload:
    """Bin valid pairs into the finest-resolution .cool (cooler cload pairix)."""
    input:
        pairs = str(RESULTS / "pairs/{sample}.valid.pairs.gz"),
        index = str(RESULTS / "pairs/{sample}.valid.pairs.gz.px2"),
        cs    = CHROMSIZES
    output:
        cool = temp(str(RESULTS / "matrix/cool/{sample}.base.cool"))
    params:
        base = BASE_RES,
        asm  = GENOME
    log:
        str(LOGS / "matrix/{sample}.cload.log")
    container:
        IMG_COOLTOOLS
    threads: 8
    resources:
        mem_mb = 32000, runtime = 480
    shell:
        r"""
        mkdir -p $(dirname {output.cool})
        cooler cload pairix -p {threads} --assembly {params.asm} \
            {input.cs}:{params.base} {input.pairs} {output.cool} > {log} 2>&1
        """


rule cooler_zoomify:
    """Coarsen to the resolution pyramid and ICE-balance each level -> .mcool."""
    input:
        cool = str(RESULTS / "matrix/cool/{sample}.base.cool")
    output:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool")
    params:
        res       = ",".join(str(r) for r in RES_LIST),
        max_iters = config["matrix"].get("balance_max_iters", 1000)
    log:
        str(LOGS / "matrix/{sample}.zoomify.log")
    container:
        IMG_COOLTOOLS
    threads: 8
    resources:
        mem_mb = 48000, runtime = 720
    shell:
        r"""
        cooler zoomify -n {threads} -r {params.res} \
            --balance --balance-args "--max-iters {params.max_iters} --nproc {threads}" \
            -o {output.mcool} {input.cool} > {log} 2>&1
        """


rule juicer_pairs:
    """
    Convert the validated .pairs to Juicer "short" format for `juicer_tools pre`.

    WHY NOT feed the .pairs directly: juicer_tools 2.20 crashes parsing the 4DN
    pairs format — `parseDCICFormat` throws ArrayIndexOutOfBoundsException
    (Index 3 / length 2), and this is INDEPENDENT of the column count (reproduced
    with both the full 10-col file and a stripped 7-col file — identical crash).
    Juicer's "short" format uses a separate, robust parser, sidestepping it.

    Short format (aidenlab/juicer Pre wiki):
        <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>
    space-separated, NO header. strand: 0=forward, non-zero=reverse (pairs use
    +/-). frag has no meaning without a restriction map, so use the required
    DISTINCT dummies frag1=0, frag2=1. Lossless for a binned .hic: the pairs are
    already MAPQ>=30 (pairtools --min-mapq 30) + UU/UR/RU (pairs_select) filtered,
    so only chrom/pos/strand matter. Run in IMG_4DN (has zcat/awk/bgzip).
    """
    input:
        pairs = str(RESULTS / "pairs/{sample}.valid.pairs.gz")
    output:
        short = temp(str(RESULTS / "pairs/{sample}.juicer.short.gz"))
    log:
        str(LOGS / "matrix/{sample}.juicer_pairs.log")
    container:
        IMG_4DN
    threads: 1
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        r"""
        ( zcat {input.pairs} \
            | awk 'BEGIN{{OFS=" "}} !/^#/ {{
                     s1 = ($6=="+") ? 0 : 1
                     s2 = ($7=="+") ? 0 : 1
                     print s1, $2, $3, 0, s2, $4, $5, 1
                   }}' \
            | bgzip > {output.short}
        ) > {log} 2>&1
        """


rule juicer_hic:
    """
    Build the Juicer .hic from the juicer "short"-format contacts (see juicer_pairs).
    juicer_tools pre computes VC, VC_SQRT, KR and SCALE normalisation vectors,
    which Arrowhead / HiCCUPS then read via `-k SCALE`.
    """
    input:
        pairs = str(RESULTS / "pairs/{sample}.juicer.short.gz"),
        cs    = CHROMSIZES,
        jar   = JUICER_JAR
    output:
        hic = str(RESULTS / "matrix/hic/{sample}.hic")
    params:
        res = ",".join(str(r) for r in HIC_RES),
        mem = config["matrix"].get("hic_java_mem", "64g")
    log:
        str(LOGS / "matrix/{sample}.juicer_pre.log")
    container:
        IMG_JAVA
    threads: 8
    resources:
        mem_mb = 80000, runtime = 1440
    shell:
        r"""
        mkdir -p $(dirname {output.hic})
        java -Xmx{params.mem} -jar {input.jar} pre \
            -r {params.res} \
            {input.pairs} {output.hic} {input.cs} > {log} 2>&1
        """
