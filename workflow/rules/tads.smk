# =============================================================================
# tads.smk — TAD / boundary calling
#
#   PRIMARY     : cooltools insulation (Crane et al. 2015 insulation score) on the
#                 ICE-balanced .mcool. Simple, reproducible, resolution-robust ->
#                 the de-facto community standard for boundaries.
#   ORTHOGONAL  : Juicer Arrowhead on the .hic. Detects (nested) contact domains
#                 via the corner score; good cross-check, widely cited (Rao 2014).
#   THIRD       : HiCExplorer hicFindTADs (Ramirez et al. 2018) on the .mcool.
#                 Insulation-based with multiple-testing-controlled boundary calls
#                 (FDR over diamond windows). Adds a packaged, FDR-aware second
#                 opinion to the cooltools score. OPT-IN via `tads.hicexplorer:true`.
#
#   boundary_dynamics.tsv summarises which boundaries appear at which timepoints.
# =============================================================================

rule insulation:
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool"),
        view  = VIEW
    output:
        tsv = str(RESULTS / "tads/insulation/{sample}.insulation.tsv")
    params:
        res     = config["tads"]["insulation_resolution"],
        windows = " ".join(str(w) for w in config["tads"]["insulation_windows"]),
        thr     = config["tads"]["insulation_threshold"],
        idiag   = config["tads"]["ignore_diags"]
    log:
        str(LOGS / "tads/{sample}.insulation.log")
    container:
        IMG_COOLTOOLS
    threads: 8
    resources:
        mem_mb = 32000, runtime = 240
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        cooltools insulation \
            {input.mcool}::resolutions/{params.res} {params.windows} \
            --view {input.view} \
            --threshold {params.thr} \
            --ignore-diags {params.idiag} \
            --clr-weight-name weight \
            -p {threads} \
            -o {output.tsv} > {log} 2>&1
        """


rule insulation_to_tads:
    input:
        tsv = str(RESULTS / "tads/insulation/{sample}.insulation.tsv")
    output:
        tads  = str(RESULTS / "tads/insulation/{sample}.tads.bed"),
        bound = str(RESULTS / "tads/insulation/{sample}.boundaries.bed")
    params:
        window = config["tads"]["insulation_primary_window"],
        chroms = ",".join(ANALYSIS_CHROMS)
    log:
        str(LOGS / "tads/{sample}.ins2tads.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30
    shell:
        r"""
        python workflow/scripts/insulation_to_tads.py \
            --insulation {input.tsv} \
            --window-bp {params.window} \
            --keep-chroms {params.chroms} \
            --out-tads {output.tads} \
            --out-boundaries {output.bound} > {log} 2>&1
        """


rule arrowhead:
    input:
        hic = str(RESULTS / "matrix/hic/{sample}.hic"),
        jar = JUICER_JAR
    output:
        domains = str(RESULTS / "tads/arrowhead/{sample}.domains.bed")
    params:
        res  = config["tads"]["arrowhead_resolution"],
        norm = config["tads"]["arrowhead_norm"],
        mem  = config["tads"]["arrowhead_java_mem"],
        outdir = lambda wc: str(RESULTS / "tads/arrowhead" / wc.sample)
    log:
        str(LOGS / "tads/{sample}.arrowhead.log")
    container:
        IMG_JAVA
    threads: 4
    resources:
        mem_mb = 72000, runtime = 480
    shell:
        r"""
        mkdir -p {params.outdir}
        java -Xmx{params.mem} -jar {input.jar} arrowhead \
            --threads {threads} -r {params.res} -k {params.norm} \
            {input.hic} {params.outdir} > {log} 2>&1
        # Arrowhead writes <outdir>/<res>_blocks.bedpe with a header row
        #   chr1  x1  x2  chr2  y1  y2  color  corner_score ...
        # Drop comment lines AND the header (identified by x1 not being numeric);
        # keep every real domain INCLUDING those on chr1.  [fix: the old filter
        # `$1!="chr1"` deleted all chromosome-1 domains, not just the header.]
        awk 'BEGIN{{OFS="\t"}} !/^#/ && $2 ~ /^[0-9]+$/ {{print $1,$2,$3}}' \
            {params.outdir}/{params.res}_blocks.bedpe \
            | sort -k1,1 -k2,2n | uniq > {output.domains}
        """


def _arrowhead_domains_for_all():
    return expand(str(RESULTS/"tads/arrowhead/{sample}.domains.bed"), sample=SAMPLES)


# ----------------------------------------------------------------------------
# Third caller: HiCExplorer hicFindTADs (Ramirez et al., 2018, Nat Commun).
# Reads the .mcool at the configured resolution via the cooler URI syntax
# (`<mcool>::/resolutions/<bp>`). Outputs (under tads/hicexplorer/):
#   {sample}.domains.bed     -- TAD calls (the canonical output)
#   {sample}.boundaries.bed  -- boundary regions
#   {sample}.score.bedgraph  -- TAD-separation score track
# Depth & step are derived from `hicexplorer_resolution` × the *_mult knobs
# in config; defaults (3*res, 10*res, 1*res) follow the HiCExplorer manual.
# `--correctForMultipleTesting fdr` is what differentiates this from the bare
# cooltools insulation: HiCExplorer FDR-controls boundaries across the swept
# diamond windows, which is robust on moderately-sparse data.
# Off by default — toggle via `tads.hicexplorer: true` in the config.
# ----------------------------------------------------------------------------
rule hicfindtads:
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool")
    output:
        domains    = str(RESULTS / "tads/hicexplorer/{sample}.domains.bed"),
        boundaries = str(RESULTS / "tads/hicexplorer/{sample}.boundaries.bed"),
        score      = str(RESULTS / "tads/hicexplorer/{sample}.score.bedgraph")
    params:
        res       = config["tads"]["hicexplorer_resolution"],
        min_depth = config["tads"]["hicexplorer_resolution"] * config["tads"]["hicexplorer_min_depth_mult"],
        max_depth = config["tads"]["hicexplorer_resolution"] * config["tads"]["hicexplorer_max_depth_mult"],
        step      = config["tads"]["hicexplorer_resolution"] * config["tads"]["hicexplorer_step_mult"],
        fdr       = config["tads"]["hicexplorer_fdr"],
        delta     = config["tads"]["hicexplorer_delta"],
        outdir    = str(RESULTS / "tads/hicexplorer"),
        prefix    = lambda wc: str(RESULTS / "tads/hicexplorer" / f"{wc.sample}.hicfind")
    log:
        str(LOGS / "tads/{sample}.hicfindtads.log")
    container:
        IMG_HICEXPLORER
    threads: 8
    resources:
        mem_mb = 32000, runtime = 480
    shell:
        r"""
        mkdir -p {params.outdir}
        hicFindTADs \
            --matrix {input.mcool}::/resolutions/{params.res} \
            --outPrefix {params.prefix} \
            --correctForMultipleTesting fdr \
            --minDepth {params.min_depth} \
            --maxDepth {params.max_depth} \
            --step {params.step} \
            --thresholdComparisons {params.fdr} \
            --delta {params.delta} \
            --numberOfProcessors {threads} > {log} 2>&1
        # hicFindTADs writes (PREFIX = {params.prefix}):
        #   <PREFIX>_domains.bed     -- TAD calls
        #   <PREFIX>_boundaries.bed  -- boundary midpoints/intervals
        #   <PREFIX>_score.bedgraph  -- per-bin TAD-separation score
        #   <PREFIX>_zscore_matrix.h5, <PREFIX>_tad_score.bm (auxiliary)
        # Normalise to the canonical output names for downstream/QC consumers.
        cp -f {params.prefix}_domains.bed    {output.domains}
        cp -f {params.prefix}_boundaries.bed {output.boundaries}
        cp -f {params.prefix}_score.bedgraph {output.score}
        """


# ----------------------------------------------------------------------------
# Cross-caller boundary agreement (per sample) + cross-sample summary.
# Compares the boundary sets from insulation / Arrowhead / hicFindTADs (the
# third one only when `tads.hicexplorer: true`); two boundaries "agree" when
# their positions sit within ±`tads.agreement_slop_bp` (default = the
# insulation_resolution). Outputs per sample:
#   tads/agreement/{sample}.agreement.tsv   -- n boundaries, pairwise Jaccards,
#                                              consensus support counts,
#                                              frac-in-consensus per caller
#   tads/agreement/{sample}.consensus.bed   -- BED5 of boundary clusters called
#                                              by >=2 callers; name = caller
#                                              list, score = caller count
# And cross-sample:
#   tads/agreement/summary.tsv              -- one row per sample (same columns)
# ----------------------------------------------------------------------------
def _agreement_hic_input(wc):
    """hicFindTADs boundaries are an input only when hicexplorer is enabled."""
    if config["tads"].get("hicexplorer", False):
        return str(RESULTS / "tads/hicexplorer" / f"{wc.sample}.boundaries.bed")
    return []


rule tad_caller_agreement:
    input:
        ins = lambda wc: str(RESULTS / "tads/insulation" / f"{wc.sample}.boundaries.bed"),
        arr = lambda wc: str(RESULTS / "tads/arrowhead"  / f"{wc.sample}.domains.bed"),
        hic = _agreement_hic_input,
    output:
        tsv       = str(RESULTS / "tads/agreement/{sample}.agreement.tsv"),
        consensus = str(RESULTS / "tads/agreement/{sample}.consensus.bed")
    params:
        slop     = config["tads"].get("agreement_slop_bp",
                                      config["tads"]["insulation_resolution"]),
        chroms   = ",".join(ANALYSIS_CHROMS),
        hic_flag = lambda wc: (
            "--hicexplorer " + str(RESULTS / "tads/hicexplorer" / f"{wc.sample}.boundaries.bed")
            if config["tads"].get("hicexplorer", False) else ""
        ),
    log:
        str(LOGS / "tads/{sample}.agreement.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python workflow/scripts/tad_caller_agreement.py \
            --sample {wildcards.sample} \
            --insulation {input.ins} \
            --arrowhead  {input.arr} \
            {params.hic_flag} \
            --slop-bp {params.slop} \
            --keep-chroms {params.chroms} \
            --out-tsv {output.tsv} \
            --out-consensus {output.consensus} > {log} 2>&1
        """


rule tad_caller_agreement_summary:
    """Cross-sample summary: concatenate the per-sample agreement TSVs."""
    input:
        per_sample = expand(str(RESULTS / "tads/agreement/{sample}.agreement.tsv"),
                            sample=SAMPLES)
    output:
        tsv = str(RESULTS / "tads/agreement/summary.tsv")
    log:
        str(LOGS / "tads/agreement_summary.log")
    threads: 1
    resources:
        mem_mb = 1000, runtime = 10
    shell:
        r"""
        # Keep header from the first per-sample TSV, append data rows from the rest.
        ( first=1
          for f in {input.per_sample}; do
              if [ "$first" = "1" ]; then
                  cat "$f"; first=0
              else
                  tail -n +2 "$f"
              fi
          done ) > {output.tsv} 2> {log}
        """


rule boundary_dynamics:
    """Union of insulation boundaries across timepoints with per-sample presence."""
    input:
        bounds = expand(str(RESULTS/"tads/insulation/{sample}.boundaries.bed"), sample=SAMPLES)
    output:
        tsv = str(RESULTS / "tads/insulation/boundary_dynamics.tsv")
    params:
        samples = ",".join(SAMPLES),
        res     = config["tads"]["insulation_resolution"]
    log:
        str(LOGS / "tads/boundary_dynamics.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30
    shell:
        r"""
        python workflow/scripts/boundary_dynamics.py \
            --beds {input.bounds} \
            --names {params.samples} \
            --merge-window {params.res} \
            --out {output.tsv} > {log} 2>&1
        """
