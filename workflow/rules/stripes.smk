################################################################################
# stripes.smk — Stripenn architectural-stripe detection.
#
# Stripenn (Yoon et al. 2022, Nat Commun) detects "architectural stripes" —
# asymmetric contact patterns extending horizontally or vertically from a single
# anchor — which are a signature of cohesin-mediated loop extrusion. They are
# distinct from focal loops (chromosight / Mustache / HiCCUPS) and from TADs.
#
# Cohesin biology specifically: stripes are interpreted as snapshots of an
# ongoing extrusion process where one anchor is "stalled" (typically at a
# CTCF site) and the cohesin ring is still extruding the loop on one side.
#
# Inputs:  the cooler-balanced .mcool file we already build.
# Outputs:
#   stripes/<sample>/result_filtered.tsv      Stripenn-filtered stripes
#   stripes/<sample>/result_unfiltered.tsv    raw stripes (pre-filter)
#
# Gated by config.stripes.enabled: true.  Stripenn is CPU-only (Python).
#
# RSchwarz 2026-06-10
################################################################################


def _stripes_cfg(key, default=None):
    return config.get("stripes", {}).get(key, default)


rule stripenn_compute:
    """
    Run Stripenn 'compute' on a sample's .mcool.  Emits a directory of per-
    chromosome results plus the headline filtered/unfiltered TSVs.
    """
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool")
    output:
        filtered   = str(RESULTS / "stripes/{sample}/result_filtered.tsv"),
        unfiltered = str(RESULTS / "stripes/{sample}/result_unfiltered.tsv"),
    params:
        res        = _stripes_cfg("resolution", 10000),
        norm       = _stripes_cfg("norm", "KR"),
        chroms     = ",".join(GREF.get("analysis_chroms", [])) or "all",
        canny_thr  = _stripes_cfg("canny", 2.0),
        max_pixel  = _stripes_cfg("max_pixel", 0.95),
        outdir     = lambda wc: str(RESULTS / "stripes" / wc.sample)
    log:
        str(LOGS / "stripes/{sample}.stripenn.log")
    container:
        IMG_STRIPENN
    threads: 4
    resources:
        mem_mb = 24000, runtime = 360
    shell:
        r"""
        mkdir -p {params.outdir}

        # Stripenn supports the cooler `.mcool` directly via the
        #   <file>::resolutions/<res> URI.
        # It writes result_filtered.tsv + result_unfiltered.tsv into --out.
        stripenn compute \
            --cool        "{input.mcool}::resolutions/{params.res}" \
            --out         "{params.outdir}/" \
            --norm        "{params.norm}" \
            --chrom       "{params.chroms}" \
            --canny       {params.canny_thr} \
            --maxpixel    {params.max_pixel} \
            --numcores    {threads} \
            > {log} 2>&1

        # Stripenn writes "result_filtered.tsv" by default; just sanity-check.
        test -s {output.filtered}
        test -s {output.unfiltered}
        """
