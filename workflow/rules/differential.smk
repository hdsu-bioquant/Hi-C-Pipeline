################################################################################
# differential.smk — replicate-aware differential Hi-C analysis.
#
# Three tools, all gated on the presence of a `condition` column in samples.tsv
# AND >= 2 conditions, each with >= 2 distinct sample names:
#
#   - HiCRep (Yang et al. 2017)             pairwise sample similarity
#                                           (stratum-adjusted correlation)
#   - multiHiCcompare (Stansfield 2018)     differential bin pairs with reps
#   - dcHiC (Chakraborty 2022)              differential compartments + TADs
#
# All three are opt-in via differential.{hicrep,multihiccompare,dchic}: true.
# They share the persistent .R_libs/ pattern (BiocManager installs on first run,
# cached for subsequent runs) — same idea as the ChIP-Seq DiffBind module.
#
# RSchwarz 2026-06-10
################################################################################

# Persistent R / Python user library lives at the pipeline root and is bind-
# mounted into the containers automatically (already covered by the SLURM
# profile's bind list).
_RLIBS = str(Path(workflow.basedir).parent / ".R_libs")


def _diff_cfg(key, default=None):
    return config.get("differential", {}).get(key, default)


# ============================================================================
# HiCRep — pairwise stratum-adjusted correlation between samples.
#
# This works on ANY pair of .mcool files, replicate or not, so we run it
# whenever the user enables it.  Output is a TSV of (sample_A, sample_B, scc).
# ============================================================================

rule hicrep_pair:
    """
    Pairwise HiCRep stratum-adjusted correlation between two .mcool files.
    Uses `hicreppy` (Python implementation).  Resolution + smoothing window are
    config-tunable.
    """
    input:
        mcool_a = str(RESULTS / "matrix/cool/{a}.mcool"),
        mcool_b = str(RESULTS / "matrix/cool/{b}.mcool"),
    output:
        tsv = str(RESULTS / "differential/hicrep/{a}_vs_{b}.tsv")
    params:
        res     = _diff_cfg("hicrep_resolution", 50000),
        h       = _diff_cfg("hicrep_smoothing_h", 5),
        max_dist = _diff_cfg("hicrep_max_distance_bp", 5_000_000),
        chroms  = ",".join(GREF.get("analysis_chroms", []))
    log:
        str(LOGS / "differential/hicrep/{a}_vs_{b}.log")
    container:
        IMG_COOLTOOLS                       # has cooler + numpy/scipy/pandas
    threads: 2
    resources:
        mem_mb = 16000, runtime = 120
    wildcard_constraints:
        a = "|".join(re.escape(s) for s in SAMPLES),
        b = "|".join(re.escape(s) for s in SAMPLES)
    shell:
        r"""
        set -e
        mkdir -p $(dirname {output.tsv})

        # Install hicreppy on first run into a persistent user-package dir.
        export PYTHONUSERBASE={_RLIBS}/pyuser
        mkdir -p $PYTHONUSERBASE
        python3 -c "import hicrep" 2>/dev/null \
            || pip install --user --no-cache-dir hicrep \
                 > {log} 2>&1 || true

        python3 - <<'PY' >> {log} 2>&1
import sys
sys.path.insert(0, "{_RLIBS}/pyuser/lib/python3.11/site-packages")
import cooler
try:
    from hicrep.utils import readMcool
    from hicrep import hicrepSCC
except Exception:
    # hicreppy 0.1.x API — alternative entry point
    from hicrep import sccByDiag as hicrepSCC
    readMcool = None

res   = {params.res}
h     = {params.h}
dBPMax = {params.max_dist}
chrom_keep = [c for c in "{params.chroms}".split(",") if c]

clrA = cooler.Cooler("{input.mcool_a}::/resolutions/" + str(res))
clrB = cooler.Cooler("{input.mcool_b}::/resolutions/" + str(res))

scc_list = hicrepSCC(clrA, clrB, h=h, dBPMax=dBPMax,
                     bDownSample=False, chromosomes=chrom_keep or None)

import numpy as np
scc_arr = np.asarray(scc_list, dtype=float)
mean_scc = float(np.nanmean(scc_arr))

with open("{output.tsv}", "w") as fh:
    fh.write("sample_a\tsample_b\tresolution\tsmoothing_h\tmax_distance_bp\tmean_scc\tper_chrom_scc\n")
    fh.write(f"{{wildcards.a}}\t{{wildcards.b}}\t{{res}}\t{{h}}\t{{dBPMax}}\t{{mean_scc:.4f}}\t{{','.join(f'{{x:.4f}}' for x in scc_arr)}}\n")

print(f"[hicrep] {{wildcards.a}} vs {{wildcards.b}}  mean_scc={{mean_scc:.4f}}")
PY
        test -s {output.tsv}
        """


rule hicrep_summary:
    """All-pairs SCC matrix — concatenate every hicrep_pair output."""
    input:
        per_pair = [str(RESULTS / "differential/hicrep" / f"{a}_vs_{b}.tsv")
                    for a, b in _combo(SAMPLES, 2)]
    output:
        summary = str(RESULTS / "differential/hicrep/summary.tsv")
    threads: 1
    resources:
        mem_mb = 1000, runtime = 10
    shell:
        r"""
        mkdir -p $(dirname {output.summary})
        head -1 {input.per_pair[0]} > {output.summary}
        for f in {input.per_pair}; do
            tail -n +2 "$f" >> {output.summary}
        done
        """


# ============================================================================
# multiHiCcompare — replicate-aware differential bin pairs.
# Requires >= 2 conditions x >= 2 replicates each.
# ============================================================================

rule multihiccompare_contrast:
    """
    Per-condition-pair differential analysis with multiHiCcompare.
    Inputs: every replicate (sample) of conditions A and B.
    Output: one TSV of differential bin pairs per contrast.
    """
    input:
        mcools_a = lambda wc: [str(RESULTS / "matrix/cool" / f"{s}.mcool")
                                for s in CONDITION_TO_SAMPLES.get(wc.condA, [])],
        mcools_b = lambda wc: [str(RESULTS / "matrix/cool" / f"{s}.mcool")
                                for s in CONDITION_TO_SAMPLES.get(wc.condB, [])],
    output:
        tsv = str(RESULTS / "differential/multihiccompare/{condA}_vs_{condB}.diff.tsv")
    params:
        res    = _diff_cfg("multihiccompare_resolution", 100000),
        fdr    = _diff_cfg("multihiccompare_fdr", 0.05),
        rlibs  = _RLIBS,
        script = str(SCRIPTS / "run_multihiccompare.R")
    log:
        str(LOGS / "differential/multihiccompare/{condA}_vs_{condB}.log")
    container:
        IMG_RDIFF
    threads: 4
    resources:
        mem_mb = 32000, runtime = 240
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        Rscript --vanilla {params.script} \
            --rlibs       {params.rlibs} \
            --mcools_a    "{input.mcools_a}" \
            --mcools_b    "{input.mcools_b}" \
            --resolution  {params.res} \
            --fdr         {params.fdr} \
            --condA       {wildcards.condA} \
            --condB       {wildcards.condB} \
            --out_tsv     {output.tsv} \
            2> {log}
        """


# ============================================================================
# dcHiC — replicate-aware differential compartments.
# Requires >= 2 conditions x >= 2 replicates each.
# ============================================================================

rule dchic_contrast:
    """
    Per-condition-pair differential compartment analysis with dcHiC.
    Inputs: every replicate (sample) .mcool for conditions A and B.
    Output: differential E1 calls + A<->B transition table.
    """
    input:
        mcools_a = lambda wc: [str(RESULTS / "matrix/cool" / f"{s}.mcool")
                                for s in CONDITION_TO_SAMPLES.get(wc.condA, [])],
        mcools_b = lambda wc: [str(RESULTS / "matrix/cool" / f"{s}.mcool")
                                for s in CONDITION_TO_SAMPLES.get(wc.condB, [])],
        chrom_sizes = CHROMSIZES,
    output:
        diff_tsv  = str(RESULTS / "differential/dchic/{condA}_vs_{condB}/dchic_diff_compartments.tsv"),
        switch_tsv = str(RESULTS / "differential/dchic/{condA}_vs_{condB}/dchic_AB_switches.tsv"),
    params:
        res    = _diff_cfg("dchic_resolution", 100000),
        rlibs  = _RLIBS,
        outdir = lambda wc: str(RESULTS / "differential/dchic" / f"{wc.condA}_vs_{wc.condB}"),
        script = str(SCRIPTS / "run_dchic.R")
    log:
        str(LOGS / "differential/dchic/{condA}_vs_{condB}.log")
    container:
        IMG_RDIFF
    threads: 4
    resources:
        mem_mb = 48000, runtime = 480
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript --vanilla {params.script} \
            --rlibs        {params.rlibs} \
            --mcools_a     "{input.mcools_a}" \
            --mcools_b     "{input.mcools_b}" \
            --condA        {wildcards.condA} \
            --condB        {wildcards.condB} \
            --resolution   {params.res} \
            --chrom_sizes  {input.chrom_sizes} \
            --out_dir      {params.outdir} \
            2> {log}
        """
