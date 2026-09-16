# =============================================================================
# loops.smk — chromatin loop / focal-peak calling (Juicer HiCCUPS)
#
#   HiCCUPS     : the original, most-validated loop caller (Rao & Huntley 2014).
#                 Local-enrichment test vs four neighbourhoods (donut, etc.).
#   APA         : Aggregate Peak Analysis -> a single QC number that confirms the
#                 called loops are, on average, real focal enrichments.
#   HiCCUPSDiff : differential loops between two timepoints.
#
# norm = SCALE (balanced); resolutions 5 kb + 10 kb (standard for mammalian loops).
# =============================================================================

rule hiccups:
    input:
        hic = str(RESULTS / "matrix/hic/{sample}.hic"),
        jar = JUICER_JAR
    output:
        merged = str(RESULTS / "loops/{sample}/merged_loops.bedpe")
    params:
        res  = ",".join(str(r) for r in config["loops"]["resolutions"]),
        norm = config["loops"]["norm"],
        mem  = config["loops"]["hiccups_java_mem"],
        outdir = lambda wc: str(RESULTS / "loops" / wc.sample / "hiccups")
    log:
        str(LOGS / "loops/{sample}.hiccups.log")
    container:
        IMG_JAVA
    threads: 8
    resources:
        mem_mb = 110000, runtime = 720
    shell:
        r"""
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        java -Xmx{params.mem} -jar {input.jar} hiccups \
            -r {params.res} -k {params.norm} \
            {input.hic} {params.outdir} > {log} 2>&1
        # juicer_tools writes the merged set as merged_loops or merged_loops.bedpe
        # depending on the build; copy whichever exists.
        if [ -s {params.outdir}/merged_loops.bedpe ]; then
            cp {params.outdir}/merged_loops.bedpe {output.merged}
        elif [ -s {params.outdir}/merged_loops ]; then
            cp {params.outdir}/merged_loops {output.merged}
        else
            echo "ERROR: HiCCUPS produced no merged_loops[.bedpe] in {params.outdir}" >&2
            echo "(usually means too few contacts for loop calling at these resolutions)" >&2
            exit 1
        fi
        """


rule apa:
    input:
        hic    = str(RESULTS / "matrix/hic/{sample}.hic"),
        merged = str(RESULTS / "loops/{sample}/merged_loops.bedpe"),
        jar    = JUICER_JAR
    output:
        done = str(RESULTS / "loops/{sample}/apa/apa.done")
    params:
        res  = config["loops"]["apa_resolution"],
        norm = config["loops"]["norm"],
        mem  = config["loops"]["apa_java_mem"],
        outdir = lambda wc: str(RESULTS / "loops" / wc.sample / "apa")
    log:
        str(LOGS / "loops/{sample}.apa.log")
    container:
        IMG_JAVA
    threads: 2
    resources:
        mem_mb = 20000, runtime = 240
    shell:
        r"""
        mkdir -p {params.outdir}
        java -Xmx{params.mem} -jar {input.jar} apa \
            -r {params.res} -k {params.norm} \
            {input.hic} {input.merged} {params.outdir} > {log} 2>&1
        touch {output.done}
        """


rule hiccupsdiff:
    input:
        hicA = lambda wc: str(RESULTS / "matrix/hic" / f"{wc.a}.hic"),
        hicB = lambda wc: str(RESULTS / "matrix/hic" / f"{wc.b}.hic"),
        loopA = lambda wc: str(RESULTS / "loops" / wc.a / "merged_loops.bedpe"),
        loopB = lambda wc: str(RESULTS / "loops" / wc.b / "merged_loops.bedpe"),
        jar  = JUICER_JAR
    output:
        done = str(RESULTS / "loops/diff/{a}_vs_{b}/differentialLoops.done")
    params:
        res  = config["loops"]["diff_resolution"],
        norm = config["loops"]["norm"],
        mem  = config["loops"]["diff_java_mem"],
        # HiCCUPSDiff wants the hiccups OUTPUT DIRS that hold the merged loop lists
        dirA = lambda wc: str(RESULTS / "loops" / wc.a / "hiccups"),
        dirB = lambda wc: str(RESULTS / "loops" / wc.b / "hiccups"),
        outdir = lambda wc: str(RESULTS / "loops/diff" / f"{wc.a}_vs_{wc.b}")
    wildcard_constraints:
        a = "|".join(SAMPLES),
        b = "|".join(SAMPLES)
    log:
        str(LOGS / "loops/{a}_vs_{b}.hiccupsdiff.log")
    container:
        IMG_JAVA
    threads: 4
    resources:
        mem_mb = 40000, runtime = 480
    shell:
        r"""
        mkdir -p {params.outdir}
        java -Xmx{params.mem} -jar {input.jar} hiccupsdiff \
            -r {params.res} -k {params.norm} \
            {input.hicA} {input.hicB} {params.dirA} {params.dirB} {params.outdir} \
            > {log} 2>&1
        touch {output.done}
        """


rule chromosight_loops:
    """
    OPTIONAL cooler-native loop caller (Matthey-Doret et al. 2020). Often more
    sensitive than HiCCUPS on moderate-depth libraries, runs on the .mcool with no
    GPU/Java. Enable with config.loops.chromosight: true. Complements (does not
    replace) HiCCUPS, which stays the most-validated reference caller.
    """
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool")
    output:
        tsv = str(RESULTS / "loops/chromosight/{sample}.chromosight.tsv")
    params:
        res    = config["loops"].get("chromosight_resolution", 10000),
        prefix = lambda wc: str(RESULTS / "loops/chromosight" / f"{wc.sample}.chromosight")
    log:
        str(LOGS / "loops/{sample}.chromosight.log")
    container:
        IMG_CHROMOSIGHT
    threads: 4
    resources:
        mem_mb = 16000, runtime = 240
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        chromosight detect --pattern loops --threads {threads} \
            {input.mcool}::resolutions/{params.res} {params.prefix} > {log} 2>&1
        # chromosight writes <prefix>.tsv
        test -s {output.tsv}
        """


# =============================================================================
# RSchwarz 2026-06-10 — Mustache loop caller + cross-caller loop agreement.
# Mustache (Roayaei Ardakany et al. 2020) is a CPU-only, cooler-native loop
# caller that uses scale-space pixel detection. It often catches loop
# populations that chromosight's template-matching misses, so running both +
# building a consensus mirrors the TAD-caller agreement pattern.
# =============================================================================

rule mustache_loops:
    """
    Mustache loop caller — scale-space pixel detection on a balanced .mcool.
    Output BEDPE has columns: chrom1, start1, end1, chrom2, start2, end2, fdr,
    detection_scale.  CPU-only.  Enable with config.loops.mustache: true.
    """
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool")
    output:
        tsv = str(RESULTS / "loops/mustache/{sample}.mustache.tsv")
    params:
        res     = config["loops"].get("mustache_resolution", 10000),
        pt      = config["loops"].get("mustache_pthreshold", 0.1),
        st      = config["loops"].get("mustache_sthreshold", 0.88),
        norm    = config["loops"].get("mustache_norm", "weight"),
        chroms  = ",".join(GREF.get("analysis_chroms", []))
    log:
        str(LOGS / "loops/{sample}.mustache.log")
    container:
        IMG_MUSTACHE
    threads: 4
    resources:
        mem_mb = 24000, runtime = 360
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        # Mustache CLI: cooler file + resolution + p-value/sparsity thresholds.
        # --norm passes through to cooler's stored weight column (ICE = 'weight').
        mustache \
            -f {input.mcool} \
            -r {params.res} \
            -pt {params.pt} \
            -st {params.st} \
            --norm {params.norm} \
            -p {threads} \
            -o {output.tsv} \
            > {log} 2>&1
        test -s {output.tsv}
        """


rule loop_caller_agreement:
    """
    Per-sample loop-caller cross-agreement.  Mirrors tad_caller_agreement:
    extract anchor pairs from each caller, compute slop-tolerant Jaccards in
    2-D (both anchors must overlap), single-linkage cluster across callers
    within the slop, emit consensus BEDPE of loops called by >= 2 callers.

    Always runs when chromosight or mustache (or both) is enabled.  HiCCUPS
    is included automatically if loops.hiccups is true.
    """
    input:
        chromosight = lambda wc: (str(RESULTS / "loops/chromosight" / f"{wc.sample}.chromosight.tsv")
                                  if config["loops"].get("chromosight", False) else []),
        mustache    = lambda wc: (str(RESULTS / "loops/mustache" / f"{wc.sample}.mustache.tsv")
                                  if config["loops"].get("mustache", False) else []),
        hiccups     = lambda wc: (str(RESULTS / "loops" / wc.sample / "merged_loops.bedpe")
                                  if config["loops"].get("hiccups", False) else [])
    output:
        agreement  = str(RESULTS / "loops/agreement/{sample}.agreement.tsv"),
        consensus  = str(RESULTS / "loops/agreement/{sample}.consensus.bedpe")
    params:
        slop_bp = config["loops"].get("agreement_slop_bp", 10000),
        script  = str(SCRIPTS / "loop_caller_agreement.py")
    log:
        str(LOGS / "loops/{sample}.loop_agreement.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p $(dirname {output.agreement})
        python3 {params.script} \
            --sample      {wildcards.sample} \
            --slop_bp     {params.slop_bp} \
            --chromosight "{input.chromosight}" \
            --mustache    "{input.mustache}" \
            --hiccups     "{input.hiccups}" \
            --out_tsv     {output.agreement} \
            --out_bedpe   {output.consensus} \
            2> {log}
        """


rule loop_caller_agreement_summary:
    """
    Cross-sample summary table — concatenate per-sample agreement rows.
    """
    input:
        per_sample = expand(str(RESULTS / "loops/agreement/{sample}.agreement.tsv"),
                            sample=SAMPLES)
    output:
        summary = str(RESULTS / "loops/agreement/summary.tsv")
    threads: 1
    resources:
        mem_mb = 1000, runtime = 10
    shell:
        r"""
        mkdir -p $(dirname {output.summary})
        # First file's header; then body of every file
        head -1 {input.per_sample[0]} > {output.summary}
        for f in {input.per_sample}; do
            tail -n +2 "$f" >> {output.summary}
        done
        """
