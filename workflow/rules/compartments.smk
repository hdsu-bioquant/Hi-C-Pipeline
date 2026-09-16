# =============================================================================
# compartments.smk — A/B compartments (cooltools eigs-cis)
#
#   eigs-cis    : cis eigendecomposition of the observed/expected matrix.
#                 E1 (first eigenvector) = the compartment profile.
#   GC phasing  : E1 sign is arbitrary; we orient it with GC content so that
#                 E1 > 0 == A (gene-rich/active). This is the accepted, genome-
#                 agnostic way to make compartments comparable across samples.
#   saddle      : compartmentalisation STRENGTH (how cleanly A-A/B-B segregate
#                 from A-B), a single comparable number per timepoint.
#   switching   : per-bin A<->B sign changes across the timecourse.
# =============================================================================

rule compartments_eigs:
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool"),
        gc    = str(RESULTS / "genome" / f"gc_{config['compartments']['resolution']}.tsv"),
        view  = VIEW
    output:
        vecs = str(RESULTS / "compartments/{sample}.cis.vecs.tsv")
    params:
        res    = config["compartments"]["resolution"],
        prefix = lambda wc: str(RESULTS / "compartments" / wc.sample)
    log:
        str(LOGS / "compartments/{sample}.eigs.log")
    container:
        IMG_COOLTOOLS
    threads: 4
    resources:
        mem_mb = 32000, runtime = 120
    shell:
        r"""
        mkdir -p $(dirname {output.vecs})
        cooltools eigs-cis \
            --phasing-track {input.gc} \
            --view {input.view} \
            --n-eigs 3 \
            -o {params.prefix} \
            {input.mcool}::resolutions/{params.res} > {log} 2>&1
        """


# The E1 -> bigWig conversion is split in two so each step runs in an image that
# actually contains its tool: pandas extraction (cooltools image) then
# bedGraphToBigWig (ucsc image). No single biocontainer has both.
rule e1_bedgraph:
    input:
        vecs = str(RESULTS / "compartments/{sample}.cis.vecs.tsv")
    output:
        bg = temp(str(RESULTS / "compartments/{sample}.E1.bedGraph"))
    log:
        str(LOGS / "compartments/{sample}.e1bg.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        # vecs.tsv columns: chrom start end ... E1 ... ; take chrom/start/end/E1, drop NaN, sort.
        python - {input.vecs} > {output.bg} 2> {log} <<'PY'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t")
df = df[["chrom","start","end","E1"]].dropna()
df = df.sort_values(["chrom","start"])
df.to_csv(sys.stdout, sep="\t", header=False, index=False)
PY
        """


rule e1_bigwig:
    input:
        bg = str(RESULTS / "compartments/{sample}.E1.bedGraph"),
        cs = CHROMSIZES
    output:
        bw = str(RESULTS / "compartments/{sample}.E1.bw")
    log:
        str(LOGS / "compartments/{sample}.e1bw.log")
    container:
        IMG_UCSC
    threads: 1
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        bedGraphToBigWig {input.bg} {input.cs} {output.bw} > {log} 2>&1
        """


rule compartment_saddle:
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool"),
        vecs  = str(RESULTS / "compartments/{sample}.cis.vecs.tsv"),
        view  = VIEW
    output:
        txt = str(RESULTS / "compartments/saddle/{sample}.saddle.txt")
    params:
        res     = config["compartments"]["resolution"],
        nbins   = config["compartments"]["n_quantile_bins"],
        extreme = config["compartments"].get("saddle_extreme_fraction", 0.2),
        prefix  = lambda wc: str(RESULTS / "compartments/saddle" / wc.sample)
    log:
        str(LOGS / "compartments/{sample}.saddle.log")
    container:
        IMG_COOLTOOLS
    threads: 4
    resources:
        mem_mb = 32000, runtime = 120
    shell:
        r"""
        mkdir -p $(dirname {output.txt})
        python workflow/scripts/compartment_saddle.py \
            --mcool {input.mcool} --resolution {params.res} \
            --eigs {input.vecs} --n-bins {params.nbins} \
            --extreme-fraction {params.extreme} --view {input.view} \
            --out-prefix {params.prefix} --out-strength {output.txt} > {log} 2>&1
        """


rule compartment_switching:
    input:
        vecs = expand(str(RESULTS/"compartments/{sample}.cis.vecs.tsv"), sample=SAMPLES)
    output:
        tsv = str(RESULTS / "compartments/compartment_switching.tsv")
    params:
        samples = ",".join(SAMPLES)
    log:
        str(LOGS / "compartments/switching.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30
    shell:
        r"""
        python workflow/scripts/compartment_switching.py \
            --vecs {input.vecs} --names {params.samples} --out {output.tsv} > {log} 2>&1
        """
