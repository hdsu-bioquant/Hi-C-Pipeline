# =============================================================================
# common.smk — config parsing, sample sheet, path helpers, target builders
#
# Adds (vs the original):
#   * fail-fast validation of genome files, the Juicer jar and every FASTQ
#   * multi-lane / multi-run support: several rows with the same `sample` are
#     aligned independently and merged at the .pairs stage (pairtools merge)
#   * an "analysis view" (chr1..22, chrX) so compartment/insulation analysis
#     ignores chrM, chrY, unplaced scaffolds and the EBV decoy
# =============================================================================
import os
import re
import sys
import pandas as pd
from pathlib import Path

# ---- Paths -----------------------------------------------------------------
RESULTS = Path(config["datapath"])
LOGS    = RESULTS / "logs"
SCRIPTS = Path(workflow.basedir) / "scripts"

GENOME      = config["genome"]
GREF        = config["genomes"][GENOME]
FASTA       = GREF["fasta"]
CHROMSIZES  = GREF["chrom_sizes"]
JUICER_ID   = GREF.get("juicer_genome_id", GENOME)
ANALYSIS_CHROMS = GREF["analysis_chroms"]
RESTRICTION_ENZYME = GREF.get("restriction_enzyme", "unspecified")

JUICER_JAR = config["juicer_tools_jar"]

# ---- Container images (Apptainer/Singularity) ------------------------------
# This pipeline runs entirely from container images (profile:
# software-deployment-method: [apptainer]); there is NO conda fallback at runtime.
# Tags below were verified to exist on quay.io / Docker Hub (2026-05).
#   * the upstream streaming step `bwa mem | pairtools parse | sort` needs bwa +
#     pairtools + samtools + pairix in ONE image -> the canonical 4D Nucleome
#     processing image (duplexa/4dn-hic) provides exactly that.
#   * cooler + cooltools (+ numpy/scipy/pandas/matplotlib) come from the cooltools
#     biocontainer (cooltools depends on, and ships, cooler).
IMG_4DN        = "docker://duplexa/4dn-hic:v43"
IMG_BWA        = "docker://quay.io/biocontainers/bwa:0.7.18--h577a1d6_2"
IMG_COOLTOOLS  = "docker://quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3"
IMG_JAVA       = "docker://eclipse-temurin:17"
IMG_FASTQC     = "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
IMG_FASTP      = "docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5"
IMG_MULTIQC    = "docker://ewels/multiqc:latest"
IMG_UCSC       = "docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:455--h2a80c09_1"
IMG_CHROMOSIGHT = "docker://quay.io/biocontainers/chromosight:1.6.3--pyhdfd78af_0"
# RSchwarz 2026-06-10 — new optional modules
IMG_MUSTACHE    = "docker://quay.io/biocontainers/mustache-hic:1.3.3--pyhdfd78af_0"
IMG_STRIPENN    = "docker://quay.io/biocontainers/stripenn:1.1.65.20--pyhdfd78af_0"
# Persistent R-library pattern for replicate-aware differential analysis
# (multiHiCcompare, dcHiC) — uses the same .R_libs/ approach as the ChIP-Seq
# pipeline's ChIPseeker/DiffBind rules.
IMG_RDIFF       = "docker://bioconductor/bioconductor_docker:RELEASE_3_19"
IMG_HICEXPLORER = "docker://quay.io/biocontainers/hicexplorer:3.7.6--pyhdfd78af_0"

# bwa index lives next to the FASTA (built by genome_prep)
BWA_PREFIX = FASTA  # bwa uses <fasta>.{amb,ann,bwt,pac,sa}

# Analysis view (BED4: chrom start end name) of the chromosomes we actually call
# features on. Built by genome_prep:analysis_view and handed to cooltools --view.
VIEW = str(RESULTS / "genome" / "analysis_view.bed")

# ---- Sample sheet ----------------------------------------------------------
# Columns (TSV): sample  timepoint  r1  r2  [lane]
#   * `sample`    biological library name (one matrix per sample)
#   * `timepoint` grouping label (free text; used for reporting/diff)
#   * `r1` / `r2` absolute paths to the paired FASTQs
#   * `lane`      OPTIONAL. If a library was sequenced over several lanes/runs,
#                 add one row per lane with the same `sample` and distinct `lane`.
#                 Lanes are aligned separately and merged before dedup.
_samples = pd.read_csv(config["samples_tsv"], sep="\t", dtype=str).fillna("")
_samples = _samples[_samples["sample"].str.strip() != ""]
if _samples.empty:
    raise SystemExit(f"[ERROR] no samples found in {config['samples_tsv']}")

SAMPLES = list(dict.fromkeys(_samples["sample"].tolist()))   # unique, order-preserving

def _sanitize(x, fallback):
    x = str(x).strip()
    x = re.sub(r"[^A-Za-z0-9._-]+", "-", x)
    return x or fallback

# SAMPLE2UNITS[sample] = [{"unit","r1","r2"}, ...]   (one entry per lane/run)
SAMPLE2UNITS = {}
SAMPLE2META  = {}
for s in SAMPLES:
    rows = _samples[_samples["sample"] == s].reset_index(drop=True)
    SAMPLE2META[s] = {"timepoint": rows.iloc[0].get("timepoint", s) or s}
    units = []
    seen = set()
    for i, row in rows.iterrows():
        unit = _sanitize(row["lane"], f"u{i+1}") if "lane" in rows.columns and str(row["lane"]).strip() else f"u{i+1}"
        if unit in seen:                       # guard against duplicate lane labels
            unit = f"{unit}-{i+1}"
        seen.add(unit)
        units.append({"unit": unit, "r1": row["r1"], "r2": row["r2"]})
    SAMPLE2UNITS[s] = units

# flat (sample, unit) list, used to expand per-lane targets (FastQC, trimming)
UNITS = [(s, u["unit"]) for s in SAMPLES for u in SAMPLE2UNITS[s]]

# ---------------------------------------------------------------------------
# Replicate / condition detection (drives the differential.smk rules).
#
# samples.tsv MAY include two extra columns:
#   condition  — free-text label grouping samples that should be compared
#                statistically (e.g. "ctrl", "treat", or a timepoint).
#   replicate  — free-text label distinguishing biological replicates within
#                a condition (e.g. "rep1", "rep2").  Used only to confirm a
#                condition has >=2 distinct replicates.
#
# If either column is absent, the replicate-aware rules silently skip
# (DIFFHIC_ELIGIBLE = False; CONDITION_TO_SAMPLES = {}).
# ---------------------------------------------------------------------------
from collections import defaultdict as _dd
from itertools import combinations as _combo

HAS_CONDITION_COL = "condition" in _samples.columns
HAS_REPLICATE_COL = "replicate" in _samples.columns

SAMPLE_TO_CONDITION = {}
CONDITION_TO_SAMPLES = _dd(list)
if HAS_CONDITION_COL:
    for s in SAMPLES:
        rows = _samples[_samples["sample"] == s].reset_index(drop=True)
        cond = str(rows.iloc[0].get("condition", "")).strip()
        if cond:
            SAMPLE_TO_CONDITION[s] = cond
            CONDITION_TO_SAMPLES[cond].append(s)

# A condition is "replicate-eligible" if it has >=2 distinct samples.
ELIGIBLE_CONDITIONS = sorted([c for c, samps in CONDITION_TO_SAMPLES.items()
                              if len(set(samps)) >= 2])
# Differential analysis needs at least 2 conditions, each replicate-eligible.
DIFFHIC_ELIGIBLE = len(ELIGIBLE_CONDITIONS) >= 2
DIFFHIC_PAIRS = list(_combo(ELIGIBLE_CONDITIONS, 2)) if DIFFHIC_ELIGIBLE else []

def timepoint(sample): return SAMPLE2META[sample]["timepoint"]

def unit_fastq(sample, unit, read):
    for u in SAMPLE2UNITS[sample]:
        if u["unit"] == unit:
            return u["r1"] if read == "r1" else u["r2"]
    raise KeyError(f"unit {unit} not found for sample {sample}")

def units_sorted_pairs(wc):
    return [str(RESULTS / "pairs" / wc.sample / f"{u['unit']}.sorted.pairs.gz")
            for u in SAMPLE2UNITS[wc.sample]]

# ---- Wildcard constraints (keep paths unambiguous with multi-lane) ----------
wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLES),
    unit   = r"[^/]+",
    r      = r"R[12]"

# ---- Resolutions -----------------------------------------------------------
RES_LIST      = config["matrix"]["resolutions"]
BASE_RES      = config["matrix"]["base_resolution"]
HIC_RES       = config["matrix"]["hic_resolutions"]

def run(step): return bool(config["runrule"].get(step, False))
def qc_on(step, default=True): return bool(config.get("qc", {}).get(step, default))

# ---- FASTQ used for alignment (raw, or fastp-trimmed if trim enabled) -------
def fq_for_align(wc):
    if run("trim"):
        return {
            "r1": str(RESULTS / "trim" / f"{wc.sample}.{wc.unit}_R1.trimmed.fastq.gz"),
            "r2": str(RESULTS / "trim" / f"{wc.sample}.{wc.unit}_R2.trimmed.fastq.gz"),
        }
    return {"r1": unit_fastq(wc.sample, wc.unit, "r1"),
            "r2": unit_fastq(wc.sample, wc.unit, "r2")}

# ---- Fail-fast input validation -------------------------------------------
def _validate_inputs():
    """Check that everything the run needs exists, with one clear error listing
    all problems, instead of failing deep inside an alignment job hours later."""
    if not config.get("validate_inputs", True):
        return
    problems = []

    # genome + chrom.sizes
    if not os.path.exists(FASTA):
        problems.append(f"genome FASTA not found: {FASTA}")
    if not os.path.exists(CHROMSIZES):
        problems.append(f"chrom.sizes not found: {CHROMSIZES}")

    # Juicer jar only required if a Juicer-dependent module is on
    needs_juicer = run("hic") or run("loops") or run("tads")
    if needs_juicer and not os.path.exists(JUICER_JAR):
        problems.append(f"juicer_tools jar not found (needed for hic/loops/tads): {JUICER_JAR}")

    # every FASTQ in the sheet
    for s in SAMPLES:
        for u in SAMPLE2UNITS[s]:
            for read in ("r1", "r2"):
                p = u[read]
                if not p:
                    problems.append(f"{s} [{u['unit']}]: empty {read} path in sample sheet")
                elif not os.path.exists(p):
                    problems.append(f"{s} [{u['unit']}]: {read} FASTQ not found: {p}")

    # analysis_chroms must exist in chrom.sizes
    if os.path.exists(CHROMSIZES):
        have = set()
        with open(CHROMSIZES) as fh:
            for line in fh:
                if line.strip():
                    have.add(line.split()[0])
        missing = [c for c in ANALYSIS_CHROMS if c not in have]
        if missing:
            problems.append(f"analysis_chroms not in chrom.sizes: {missing}")

    if problems:
        msg = "\n  - ".join(problems)
        raise SystemExit(
            f"[ERROR] input validation failed ({len(problems)} problem(s)):\n  - {msg}\n"
            f"Fix config/samples.tsv and config/config.yaml, or set validate_inputs: false to skip."
        )

_validate_inputs()

# ---- Target builder (what `rule all` asks for) -----------------------------
def targets():
    t = []
    if run("fastqc"):
        t += [str(RESULTS/f"qc/fastqc/{s}.{u}_{r}_fastqc.zip")
              for (s, u) in UNITS for r in ["R1", "R2"]]
    if run("align_pairs"):
        t += expand(str(RESULTS/"pairs/{sample}.valid.pairs.gz"), sample=SAMPLES)
        t += expand(str(RESULTS/"pairs/{sample}.dedup.stats.txt"), sample=SAMPLES)
    if run("cooler"):
        t += expand(str(RESULTS/"matrix/cool/{sample}.mcool"), sample=SAMPLES)
    if run("hic"):
        t += expand(str(RESULTS/"matrix/hic/{sample}.hic"), sample=SAMPLES)
    if run("multiqc"):
        t += [str(RESULTS/"qc/multiqc/multiqc_report.html")]
    # Hi-C-specific QC (contact-distance decay + cross-sample library summary)
    if run("cooler") and qc_on("contact_decay"):
        t += expand(str(RESULTS/"qc/scaling/{sample}.scaling.tsv"), sample=SAMPLES)
        t += [str(RESULTS/"qc/scaling/contact_decay.png")]
    if run("align_pairs") and qc_on("library_summary"):
        t += [str(RESULTS/"qc/pairs_summary.tsv")]
    if run("tads"):
        t += expand(str(RESULTS/"tads/insulation/{sample}.tads.bed"), sample=SAMPLES)
        t += expand(str(RESULTS/"tads/arrowhead/{sample}.domains.bed"), sample=SAMPLES)
        t += [str(RESULTS/"tads/insulation/boundary_dynamics.tsv")]
        # Optional third opinion: HiCExplorer hicFindTADs (insulation-based with
        # FDR-controlled boundary calls). Off by default; flip on in config to
        # cross-check the cooltools / Arrowhead boundary sets.
        if config["tads"].get("hicexplorer", False):
            t += expand(str(RESULTS/"tads/hicexplorer/{sample}.domains.bed"), sample=SAMPLES)
        # Cross-caller boundary agreement: per-sample TSV + >=2-of-N consensus BED,
        # plus a cross-sample summary. Adapts to 2-caller (insulation+Arrowhead)
        # or 3-caller mode automatically based on `tads.hicexplorer`.
        t += expand(str(RESULTS/"tads/agreement/{sample}.agreement.tsv"), sample=SAMPLES)
        t += expand(str(RESULTS/"tads/agreement/{sample}.consensus.bed"), sample=SAMPLES)
        t += [str(RESULTS/"tads/agreement/summary.tsv")]
    if run("loops"):
        # HiCCUPS + APA + HiCCUPSDiff REQUIRE a CUDA GPU (jcuda compiles a PTX
        # kernel via nvcc at startup). On CPU-only clusters set
        # `loops.hiccups: false` in the config and use the cooler-native
        # chromosight (below) — it runs entirely on the CPU.
        if config["loops"].get("hiccups", True):
            t += expand(str(RESULTS/"loops/{sample}/merged_loops.bedpe"), sample=SAMPLES)
            t += expand(str(RESULTS/"loops/{sample}/apa/apa.done"), sample=SAMPLES)
            for pair in config["loops"].get("diff_pairs", []):
                a, b = pair.split()
                t += [str(RESULTS/f"loops/diff/{a}_vs_{b}/differentialLoops.done")]
        # optional cooler-native (CPU) loop caller; recommended on GPU-less clusters
        if config["loops"].get("chromosight", False):
            t += expand(str(RESULTS/"loops/chromosight/{sample}.chromosight.tsv"), sample=SAMPLES)
        # second cooler-native (CPU) loop caller — Mustache, scale-space pixel detection
        if config["loops"].get("mustache", False):
            t += expand(str(RESULTS/"loops/mustache/{sample}.mustache.tsv"), sample=SAMPLES)
        # cross-caller loop agreement — runs when >= 2 callers are enabled
        _n_loop_callers = sum([
            bool(config["loops"].get("chromosight", False)),
            bool(config["loops"].get("mustache",    False)),
            bool(config["loops"].get("hiccups",     False)),
        ])
        if _n_loop_callers >= 2:
            t += expand(str(RESULTS/"loops/agreement/{sample}.agreement.tsv"), sample=SAMPLES)
            t += expand(str(RESULTS/"loops/agreement/{sample}.consensus.bedpe"), sample=SAMPLES)
            t += [str(RESULTS/"loops/agreement/summary.tsv")]
    if run("compartments"):
        t += expand(str(RESULTS/"compartments/{sample}.E1.bw"), sample=SAMPLES)
        t += [str(RESULTS/"compartments/compartment_switching.tsv")]
        if config["compartments"].get("saddle", False):
            t += expand(str(RESULTS/"compartments/saddle/{sample}.saddle.txt"), sample=SAMPLES)
    # --- Stripenn — architectural stripes (cohesin biology) ----------------
    if config.get("stripes", {}).get("enabled", False):
        t += expand(str(RESULTS/"stripes/{sample}/result_filtered.tsv"), sample=SAMPLES)
    # --- Replicate-aware differential analysis ----------------------------
    # HiCRep runs on any pair of samples (no replicate requirement).
    if config.get("differential", {}).get("hicrep", False) and len(SAMPLES) >= 2:
        from itertools import combinations as __cmb
        t += [str(RESULTS/"differential/hicrep" / f"{a}_vs_{b}.tsv")
              for a, b in __cmb(SAMPLES, 2)]
        t += [str(RESULTS/"differential/hicrep/summary.tsv")]
    # multiHiCcompare + dcHiC require condition column + >=2 reps per condition
    if DIFFHIC_ELIGIBLE:
        if config.get("differential", {}).get("multihiccompare", False):
            for a, b in DIFFHIC_PAIRS:
                t += [str(RESULTS/f"differential/multihiccompare/{a}_vs_{b}.diff.tsv")]
        if config.get("differential", {}).get("dchic", False):
            for a, b in DIFFHIC_PAIRS:
                t += [str(RESULTS/f"differential/dchic/{a}_vs_{b}/dchic_diff_compartments.tsv")]
                t += [str(RESULTS/f"differential/dchic/{a}_vs_{b}/dchic_AB_switches.tsv")]
    # --- HiGlass prep -----------------------------------------------------
    if config.get("higlass", {}).get("enabled", False):
        t += [str(RESULTS/"higlass/manifest.yaml")]
        t += [str(RESULTS/"higlass/README.md")]
    return t
