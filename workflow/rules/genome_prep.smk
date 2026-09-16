# =============================================================================
# genome_prep.smk — one-time reference preparation
#   * bwa index of the genome FASTA
#   * analysis "view" BED (chr1..22, chrX) handed to cooltools --view
#   * GC-content phasing track (per compartment-resolution bin) for cooltools eigs
# =============================================================================

rule bwa_index:
    """Build the bwa index next to the FASTA (idempotent; ~1-2 h for hg38)."""
    input:
        fa = FASTA
    output:
        # bwa writes <fasta>.{amb,ann,bwt,pac,sa}
        idx = multiext(FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        str(LOGS / "genome_prep/bwa_index.log")
    container:
        IMG_BWA
    threads: 1
    resources:
        mem_mb = 8000, runtime = 240
    shell:
        "bwa index {input.fa} > {log} 2>&1"


rule analysis_view:
    """
    Whole-chromosome viewframe (BED4: chrom start end name) for the chromosomes in
    config.analysis_chroms. cooltools --view uses this so insulation, eigs-cis,
    expected and saddle ignore chrM, chrY, unplaced scaffolds and the EBV decoy
    (all present in hg38.chrom.sizes), which otherwise add noise and waste compute.
    """
    input:
        cs = CHROMSIZES
    output:
        view = VIEW
    params:
        chroms = ",".join(ANALYSIS_CHROMS)
    log:
        str(LOGS / "genome_prep/analysis_view.log")
    threads: 1
    resources:
        mem_mb = 2000, runtime = 10
    run:
        sizes = {}
        with open(input.cs) as fh:
            for line in fh:
                if line.strip():
                    name, length = line.split()[:2]
                    sizes[name] = int(length)
        os.makedirs(os.path.dirname(output.view), exist_ok=True)
        with open(output.view, "w") as out:
            for c in params.chroms.split(","):
                if c not in sizes:
                    raise SystemExit(f"[ERROR] analysis chrom {c} not in {input.cs}")
                out.write(f"{c}\t0\t{sizes[c]}\t{c}\n")


rule gc_track:
    """
    GC content per bin at the compartment resolution, restricted to the analysis
    chromosomes. Used to ORIENT the compartment eigenvector so that E1>0 == A
    (active) — the genome-agnostic, accepted way to phase A/B.
    """
    input:
        fa   = FASTA,
        cs   = CHROMSIZES,
        view = VIEW
    output:
        gc  = str(RESULTS / "genome" / f"gc_{config['compartments']['resolution']}.tsv")
    params:
        res = config["compartments"]["resolution"]
    log:
        str(LOGS / "genome_prep/gc_track.log")
    container:
        IMG_COOLTOOLS
    threads: 1
    resources:
        mem_mb = 8000, runtime = 60
    shell:
        r"""
        mkdir -p $(dirname {output.gc})
        # Restrict to the analysis chromosomes by binnifying a chrom.sizes subset
        # derived from the view (cols: chrom 0 end name -> chrom<TAB>end).
        ( awk 'BEGIN{{OFS="\t"}} {{print $1,$3}}' {input.view} > {output.gc}.chromsizes
          cooltools genome binnify {output.gc}.chromsizes {params.res} > {output.gc}.bins.bed
          cooltools genome gc {output.gc}.bins.bed {input.fa} > {output.gc}
          rm -f {output.gc}.chromsizes
        ) > {log} 2>&1
        """
