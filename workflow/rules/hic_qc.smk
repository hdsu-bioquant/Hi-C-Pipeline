# =============================================================================
# hic_qc.smk — Hi-C-specific quality control (additive; reads existing outputs)
#
#   contact_decay  : P(s), the contact frequency vs genomic distance curve — the
#                    single most diagnostic Hi-C QC. A healthy mammalian library
#                    decays ~ s^-1 over 0.1-1 Mb; a flat/short-range-only or noisy
#                    curve flags under-digestion, random ligation or low complexity.
#   contact_decay_plot : overlay all samples' P(s) on one log-log figure.
#   pairs_summary  : collate the key pairtools dedup metrics (dup rate, cis/trans,
#                    long-cis fraction, valid pairs) into one cross-sample table.
# =============================================================================

rule contact_decay:
    input:
        mcool = str(RESULTS / "matrix/cool/{sample}.mcool"),
        view  = VIEW
    output:
        tsv = str(RESULTS / "qc/scaling/{sample}.scaling.tsv")
    params:
        res   = config.get("qc", {}).get("scaling_resolution", 10000),
        idiag = config["tads"].get("ignore_diags", 2)
    log:
        str(LOGS / "qc/{sample}.scaling.log")
    container:
        IMG_COOLTOOLS
    threads: 4
    resources:
        mem_mb = 16000, runtime = 120
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python workflow/scripts/contact_decay.py \
            --mcool {input.mcool} --resolution {params.res} \
            --view {input.view} --ignore-diags {params.idiag} \
            --nproc {threads} --out {output.tsv} > {log} 2>&1
        """


rule contact_decay_plot:
    input:
        tsv = expand(str(RESULTS/"qc/scaling/{sample}.scaling.tsv"), sample=SAMPLES)
    output:
        png = str(RESULTS / "qc/scaling/contact_decay.png")
    params:
        names = ",".join(SAMPLES)
    log:
        str(LOGS / "qc/contact_decay_plot.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 4000, runtime = 30
    shell:
        r"""
        python workflow/scripts/plot_scaling.py \
            --scaling {input.tsv} --names {params.names} --out {output.png} > {log} 2>&1
        """


rule pairs_summary:
    input:
        stats = expand(str(RESULTS/"pairs/{sample}.dedup.stats.txt"), sample=SAMPLES)
    output:
        tsv = str(RESULTS / "qc/pairs_summary.tsv")
    params:
        names = ",".join(SAMPLES)
    log:
        str(LOGS / "qc/pairs_summary.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 2000, runtime = 15
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        python workflow/scripts/pairs_summary.py \
            --stats {input.stats} --names {params.names} --out {output.tsv} > {log} 2>&1
        """
