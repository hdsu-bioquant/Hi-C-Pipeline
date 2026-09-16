# Hi-C pipeline — User guide

## Executive Summary

- end-to-end: paired-end FASTQ → valid pairs → contact matrices →A/B compartments, TADs, chromatin loops + QC 
- alignment and pair processing follow the **4D Nucleome standard**
- contact matrices in both community formats, built from the same valid pairs:
  ICE-balanced multi-resolution `.mcool` (cooler / cooltools) and `.hic` (Juicer Tools, SCALE / KR)
- **3 TAD callers** (cooltools insulation, Juicer Arrowhead, HiCExplorer hicFindTADs) + cross-caller
  agreement metrics + consensus boundary set
- **CPU-only loop calling** (chromosight; optional Mustache + loop-caller consensus) -> HiCCUPS / APA /HiCCUPSDiff only for GPU clusters
- **A/B compartments** by `cooltools eigs-cis`
- every rule runs inside a pinned **Apptainer / Singularity** container; jobs are dispatched to **SLURM** through a Snakemake profile
- primary outputs: `pairs/<sample>.valid.pairs.gz`, `matrix/cool/<sample>.mcool`, `matrix/hic/<sample>.hic`

```
   raw FASTQ (R1/R2)
      │  FastQC                                          [optional: fastp adapter-only trim]
      ▼
   bwa mem -SP5M ─▶ pairtools parse ─▶ pairtools sort              
      │
      ▼
   pairtools merge ─▶ dedup ─▶ select UU/UR/RU ─▶ <sample>.valid.pairs.gz  (+ pairix index)
      │
      ├── cooler cload  ─▶ cooler zoomify --balance ─▶ <sample>.mcool   (ICE, 1 kb … 1 Mb)
      └── Juicer "short" format ─▶ juicer_tools pre ─▶ <sample>.hic     (SCALE / KR)
                                                        │
   ┌──────────────────────┬──────────────────────┬──────┴───────────────┬──────────────────────┬──────────────────────┐
   ▼                      ▼                      ▼                      ▼                      ▼                      ▼
 QC                    COMPARTMENTS           TADs                   LOOPS                  OPTIONAL               VIEWER
 P(s) contact decay    cooltools eigs-cis     cooltools insulation   chromosight            Stripenn stripes       HiGlass
 pairs summary table   GC-phased E1 → bigWig  Juicer Arrowhead       Mustache               HiCRep similarity      manifest
 MultiQC (FastQC +     saddle strength        HiCExplorer FindTADs   loop-caller agreement  multiHiCcompare /
 pairtools stats)      A↔B switching table    TAD-caller agreement   HiCCUPS + APA          differential E1
```

> The pipeline is on GitHub: **https://github.com/hdsu-bioquant/Hi-C-Pipeline**.

### How the pipeline works — step by step

1. **Read QC (FastQC).** Every FASTQ is checked for base quality, adapter content and duplication. Nothing is
   modified — it is a report only. Hi-C reads are **not trimmed** by default (see §3.2 for why).
2. **Alignment (`bwa mem -SP5M`).** The two reads of a pair are aligned *independently* to the genome. A Hi-C
   read often runs through the ligation junction, so one read can map to two places; the `-5` flag keeps the
   5′ part as the primary alignment, exactly as the 4DN consortium recommends.
3. **Pair calling (`pairtools parse` + `sort`).** Each read pair becomes one line of a `.pairs` file —
   `chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type`. The pair type says how confident the two ends are
   (`UU` both unique, `UR`/`RU` one end rescued from a multi-fragment "walk", `DD` duplicate, …). A read end
   needs MAPQ ≥ 30 to count as unique.
4. **Lane merge, de-duplication, selection.** Lanes of the same library are merged, PCR/optical duplicates are
   marked (`DD`) and only `UU`/`UR`/`RU` pairs are kept → `<sample>.valid.pairs.gz` (bgzip + pairix index).
   This file is the pipeline's **source of truth**; everything downstream is derived from it.
5. **Contact matrices.** The pairs are binned at 1 kb (`cooler cload`), coarsened to a 1 kb – 1 Mb pyramid and
   ICE-balanced (`cooler zoomify --balance`) → `.mcool`. In parallel the same pairs are converted to Juicer's
   "short" format and `juicer_tools pre` builds the `.hic` with its normalisation vectors.
6. **Matrix / library QC.** A P(s) contact-decay curve per sample (+ overlay plot), a cross-sample table of
   duplicate rate / cis fraction / valid pairs, and one MultiQC report (FastQC + pairtools statistics).
7. **A/B compartments.** `cooltools eigs-cis` at 100 kb on the analysis chromosomes; E1 is oriented with the GC
   content so that **E1 > 0 = A** in every sample. Outputs: E1 bigWig (load in IGV), a saddle-plot compartment
   strength (one number per sample) and a per-bin A↔B switching table across samples.
8. **TADs.** Three independent callers — cooltools insulation (10 kb bins, 200 kb window, Li threshold),
   Juicer Arrowhead on the `.hic` (10 kb, SCALE) and HiCExplorer hicFindTADs (FDR-controlled). Their boundaries
   are then compared per sample (slop-tolerant Jaccards, ≥ 2-of-3 consensus BED) and across samples
   (boundary-dynamics table).
9. **Loops.** chromosight on the `.mcool` at 10 kb (CPU). Optionally Mustache (CPU) plus a chromosight–Mustache
   consensus. HiCCUPS (+ APA, HiCCUPSDiff) only runs on a GPU cluster.
10. **Optional modules.** Stripenn architectural stripes, HiCRep sample similarity, replicate-aware
    differential analysis, and a HiGlass viewer manifest (§8).


## Table of contents

0. Initial Setup
1. Session Setup
2. Sample Sheet construction
3. Pipeline configuration
4. Pipeline execution
5. Outputs & Interpretation
6. TAD calling
7. Loop calling
8. Optional modules
9. Reproducibility
10. References

---

## 0. Initial Setup

You only need this section the first time you ever use the pipeline. Subsequent sessions skip to §1.

### 0.1 Requirements

| Software | Explanatory Section |
|---|---|
| **Conda / Miniconda** | §0.2 |
| **Snakemake 8 or 9** | §0.3 |
| **Apptainer or Singularity** | §0.4 |
| **Kerberos ticket (+ renewer)** | §0.5 |
| **Reference genome + chrom.sizes + Juicer Tools jar** | §0.6 |
| **Sample Sheet + FASTQ files** | §2 |

The following resources must exist on disk (they are referenced by absolute path in `config/config.yaml` and
are **not** part of the repository):

| Resource | Size | Notes |
|---|---|---|
| **Genome FASTA** (hg38 primary assembly) | ~3 GB | `genomes.hg38.fasta`. The bwa index (`.amb .ann .bwt .pac .sa`) is built **next to it** once by `rule bwa_index` (~1–2 h). |
| **chrom.sizes** (2 columns) | KB | `genomes.hg38.chrom_sizes`. Must contain every chromosome listed in `analysis_chroms`. |
| **Juicer Tools jar** (2.20.00) | ~50 MB | `juicer_tools_jar`. Needed for the `.hic`, Arrowhead, HiCCUPS/APA/HiCCUPSDiff. |
| **Container images** | ~2.6 GB | 9 images for the default run, pulled automatically on first use into `.snakemake/singularity/`. |

**The following software does *not* need to be installed:**
Docker (Apptainer reads Docker images), bwa, pairtools, pairix, samtools, cooler, cooltools, chromosight,
Mustache, HiCExplorer, Stripenn, Java, FastQC, fastp, MultiQC, bedGraphToBigWig, R/Bioconductor — these are all
inside containers (§9).

### 0.2 Miniconda Installation (skip if already available)

```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p $HOME/miniconda3
rm -f /tmp/miniconda.sh
$HOME/miniconda3/bin/conda init bash
source ~/.bashrc
conda --version
```

### 0.3 Create the Snakemake pipeline environment

```bash
conda create -y -n snakemake-9.19.0 \
    -c conda-forge -c bioconda \
    'snakemake=9.19.0' \
    snakemake-executor-plugin-slurm \
    pandas
```

Notes:
- `snakemake-executor-plugin-slurm` is **required** to dispatch jobs via `--profile profiles/slurm`. Without it,
  the SLURM profile cannot submit jobs. On a different scheduler, install the matching plugin instead.
- `pandas` is imported by `workflow/rules/common.smk` to parse the sample sheet (it is usually pulled in as a
  Snakemake dependency, but listing it explicitly is safe).
- Snakemake 8 is supported as well (`'snakemake>=8'`; the Snakefile declares `min_version("8.0")`).

Verify:
```bash
conda activate snakemake-9.19.0
snakemake --version
```

### 0.4 Check the status of Apptainer / Singularity on the Cluster

```bash
module avail system/singularity
module load system/singularity/3.10.4
singularity --version
```

If you are working on the BioQuant nodes, Apptainer does not need to be installed. The BioQuant nodes already
contain it; you only need to load the module before launching the pipeline (see §1.2). On another system make
sure `apptainer` or `singularity` is on your `PATH`.

### 0.5 Set up a persistent Kerberos renewer

The pipeline is currently maintained on `bq-storage`, a **Kerberos-authenticated NFS mount**. To read or write
anything there you need a Kerberos ticket. A ticket lasts 24 h and is renewable without a password for 7 days.
A full Hi-C run on four deep libraries takes **1–2 days**, so multi-day runs need a tiny background loop that
refreshes the ticket every 6 h.

```bash
# 1. Get a long renewable ticket (you will be prompted for the password of your BIOQUANT account)
kinit -r 7d

# 2. Verify
klist
#   "Expires"     ~24 h after Valid starting    (the ticket lifetime)
#   "renew until" ~7 days after Valid starting  (the renewable window)
```

If `renew until` equals `Expires`, the ticket isn't renewable — contact the cluster admins.

**Recommendation: start a dedicated ticket renewer in a separate tmux session:**

```bash
tmux new -s krb-keeper
# inside the krb-keeper tmux:
while sleep 21600; do
    kinit -R && echo "krb renewed $(date)" \
          || echo "krb RENEW FAILED $(date)"
done
# Detach: Ctrl-b then d
```

The credential cache is a **single file inside your home directory** (e.g. `DIR::/home/$USER/krb5cc_5067_cluster/tkt`).
One renewer covers every tmux/ssh session, every SLURM job, and the Snakemake controller. Around day 6 of a
continuous run, run `kinit -r 7d` again with your password to reset the renewable window.
(Alternative used in the profile comments: `krenew -K 60 -- snakemake …` renews the ticket for the lifetime of
one command.)

**The failure signature of a lapsed ticket:** SLURM jobs die at launch with
`slurmstepd: couldn't chdir … TaskProlog failed`, exit code `0:53`, 0 s elapsed, no SLURM log written.
If you see that, run `klist` before touching any code.

### 0.6 Reference genome, chrom.sizes and Juicer Tools (one-time)

Nothing is downloaded automatically. `config/config.yaml` is shipped as a **template**: the entries marked
`<<< EDIT >>>` (results directory, genome FASTA, chrom.sizes, Juicer jar, sample sheet) must be filled in before
the first run; `validate_inputs: true` stops the pipeline immediately with a clear list of whatever is missing.
Point the config at:

```yaml
genomes:
  hg38:
    fasta:       /path/to/your/genome/GRCh38.primary_assembly.genome.fa   # bwa index is built next to it on first run
    chrom_sizes: /path/to/your/genome/hg38.chrom.sizes                    # name <TAB> length
juicer_tools_jar: /path/to/your/tools/juicer_tools.2.20.00.jar
```

The Juicer Tools jar can be downloaded from
https://github.com/aidenlab/Juicebox/releases; it needs Java 11+ (provided by the `eclipse-temurin:17` image).

---

## 1. Session Setup

### 1.1 SSH connection and tmux session activation

If you close or interrupt the ssh connection, every process whose terminal closes is killed — including a
multi-day controller. tmux makes the controller independent from your local ssh connection.

```bash
ssh your_username@appl2.bioquant.uni-heidelberg.de
# You will be prompted for the password of your BIOQUANT account

ssh your_username@login1.bioquant.uni-heidelberg.de
# You will be prompted for the password of your BIOQUANT account

tmux new -s hic     # or `tmux attach -t hic` to re-attach to the existent tmux session
```

### 1.2 Apptainer activation

```bash
module load system/singularity/3.10.4
```

### 1.3 Snakemake environment activation

```bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate snakemake-9.19.0
```

### 1.4 Navigation to the pipeline directory

```bash
cd /net/bq-storage/ag-cherrmann/projects/38_Pipelines/Hi-C
```

> **Warning**: the directory must NOT contain `&` or any whitespace. Paths are substituted unquoted into
> container/shell commands; a bare `&` is read as a background operator.

### 1.5 Sanity check: ticket + renewer alive

```bash
klist | head -3                            # ticket valid for hours/days
pgrep -af -u "$USER" 'sleep 21600'         # renewer alive (one PID)
```

If `klist` says *Credentials cache not found* or the expiry is in the past, redo §0.5 before launching.

---

## 2. Sample Sheet construction

The raw data enter the pipeline through a **sample sheet** (`config/samples.tsv`). It is **tab-separated WITH a
header row**. One row per **library × lane**. As shipped, `config/samples.tsv` is a two-row **placeholder template**
(`/path/to/your/data/…`) — replace its rows with your own libraries. Columns:

```
sample    timepoint    r1    r2    [lane]    [condition]    [replicate]
```

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | library name → one matrix per sample; must be unique across libraries; used in every output file name |
| `timepoint` | no | free-text grouping label (defaults to `sample`); used for reporting only |
| `r1`, `r2` | yes | **absolute** paths to the paired FASTQs (`.fastq.gz` / `.fq.gz`) |
| `lane` | no | for libraries sequenced over several lanes/runs: one row per lane with the same `sample` and a distinct `lane` label. Rows without a `lane` get the unit `u1`, `u2`, … |
| `condition`, `replicate` | no | only read by the replicate-aware differential module (§8.5): samples sharing a `condition` are treated as replicates |

**Example sheet** (hNCC neural-crest time-course, hg38, single lane per library):

```
sample     timepoint  r1                                                    r2
T1_hNCC    T1         /abs/.../t1_S2_L001_R1_001.fastq.gz  /abs/.../t1_S2_L001_R2_001.fastq.gz
T2_hNCC    T2         /abs/.../t2_S3_L001_R1_001.fastq.gz  /abs/.../t2_S3_L001_R2_001.fastq.gz
T3_hNCC    T3         /abs/.../t3_S4_L001_R1_001.fastq.gz  /abs/.../t3_S4_L001_R2_001.fastq.gz
T4_hNCC    T4        /abs/.../t4_S1_L001_R1_001.fastq.gz /abs/.../t4_S1_L001_R2_001.fastq.gz
```

**Multi-lane library** (same `sample`, different `lane`):

```
sample    timepoint  r1                              r2                              lane
HCT116    ctrl       /abs/.../HCT116_L001_R1.fq.gz   /abs/.../HCT116_L001_R2.fq.gz   L001
HCT116    ctrl       /abs/.../HCT116_L002_R1.fq.gz   /abs/.../HCT116_L002_R2.fq.gz   L002
```

Lanes are aligned independently (FastQC also runs per lane), then `rule merge_lanes` (`pairtools merge`)
combines the per-lane sorted pairs **before** `pairtools dedup`, so cross-lane PCR duplicates are caught.

### Notes

- FASTQ paths must be **absolute** and exist when the controller starts (enforced by `validate_inputs: true`
  at parse time, together with the genome FASTA, chrom.sizes and the Juicer jar).
- Sample names are sanitised to `[A-Za-z0-9._-]`; duplicate `lane` labels within a sample are made unique.

---

## 3. Pipeline configuration

Everything is driven by `config/config.yaml`. Fill in the **general template** in place. 

### 3.1 Main settings

```yaml
samples_tsv: config/samples.tsv                      # sample sheet (§2)
datapath:    /path/to/your/results                   # outputs destination (absolute)
validate_inputs: true                                # fail fast on missing FASTQ / genome / jar

genome: hg38
genomes:
  hg38:
    fasta:        /path/to/your/genome/GRCh38.primary_assembly.genome.fa
    chrom_sizes:  /path/to/your/genome/hg38.chrom.sizes
    juicer_genome_id: hg38                           # label written into the .hic
    restriction_enzyme: "unspecified"                # metadata only (e.g. "MboI", "Arima (DpnII+HinfI)")
    analysis_chroms: ["chr1", …, "chr22", "chrX"]    # the "analysis view"

juicer_tools_jar: /path/to/your/tools/juicer_tools.2.20.00.jar
```

`analysis_chroms` defines the **analysis view**: compartments, insulation, saddle, P(s), the TAD/loop
agreement and hicFindTADs-based statistics are restricted to these chromosomes. The `.mcool`/`.hic` matrices
themselves are still built genome-wide. Leave out chrM, chrY (if not relevant), `chrN_*_random`, `chrUn_*` and
the EBV decoy — they add noise and compute.

### 3.2 Select steps to be run

```yaml
runrule:
  fastqc        : true     # raw read QC
  trim          : false    # KEEP OFF for Hi-C (adapter-only fastp; see below)
  align_pairs   : true     # bwa mem -SP5M | pairtools → valid pairs
  cooler        : true     # ICE-balanced multi-resolution .mcool
  hic           : true     # Juicer .hic (needed by Arrowhead / HiCCUPS)
  multiqc       : true     # FastQC + pairtools stats report
  tads          : true     # insulation + Arrowhead (+ hicFindTADs) + agreement
  loops         : true     # chromosight / Mustache (CPU) or HiCCUPS (GPU)
  compartments  : true     # eigs-cis (A/B), GC-phased, + saddle + switching

qc:
  contact_decay:   true    # P(s) curve per sample + overlay plot
  library_summary: true    # cross-sample table of dup / cis / valid pairs
  scaling_resolution: 10000
```

> **Why `trim: false`?** Hi-C reads run **through the ligation junction**, so the 5′ part of a read maps to one
> locus and the 3′ part to another. Quality trimming destroys exactly that junction-spanning signal. The correct
> handling is `bwa mem -5` + `pairtools parse --walks-policy 5unique`, both of which work on the *untrimmed*
> read (4DN: Reiff *et al.* 2022). Only for libraries with verifiable adapter contamination, `trim: true` runs
> `fastp` with quality/length filtering **disabled** (adapter-only).

### 3.3 Tunable parameters

#### Alignment + pair processing

```yaml
bwa:
  flags:   "-SP5M"                  # 4DN-canonical; do not change
  threads: 16

pairtools:
  walks_policy: "5unique"           # 4DN-canonical handling of multi-fragment reads
  max_inter_align_gap: 30
  min_mapq: 30                      # 4DN / distiller / nf-core standard
  select_expr: '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")'
  sort_threads: 8
  sort_mem: "20G"
```

| Knob | Options & guidance |
|---|---|
| `bwa.flags` | `-S` skip mate rescue, `-P` skip read pairing, `-5` 5′-most alignment of a chimeric read is primary, `-M` mark shorter split hits secondary. **Keep `-SP5M`.** |
| `min_mapq` | **`30`** *(default)*. Lower to `1` only for permissive, exploratory scans. |
| `walks_policy` | **`5unique`** *(rescues junction reads)*. `mask` ignores walks; `all` keeps every alignment of a walk (noisy). |
| `select_expr` | pair types kept after dedup. Duplicates are *marked* (`DD`) by `pairtools dedup --mark-dups` and discarded here. |
| `sort_mem` | **`"20G"`** for ~300 M-pair libraries; `"40G"` for very deep libraries. |

#### Contact matrices

```yaml
matrix:
  base_resolution:    1000
  resolutions:        [1000, 2000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
  hic_resolutions:    [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
  balance_max_iters:  1000
  hic_java_mem:       "64g"
```

| Knob | Options |
|---|---|
| `base_resolution` | **`1000`** (1 kb). `500` only for very deep libraries. |
| `resolutions` | the `.mcool` pyramid; each level is ICE-balanced (`weight` column). Loops live at 5–10 kb, TAD boundaries ~10 kb, compartments 100–250 kb. |
| `hic_resolutions` | the `.hic` levels (Juicer `pre` does not take 2 kb by default). |
| `balance_max_iters` | ICE iteration cap (typically converges in 50–200). |
| `hic_java_mem` | Java heap for `juicer_tools pre`; `32g` is enough for shallow libraries. |

### 3.4 TAD-caller settings

```yaml
tads:
  insulation_resolution: 10000
  insulation_windows:    [100000, 200000, 500000]
  insulation_primary_window: 200000
  insulation_threshold:  "Li"            # Li | Otsu
  ignore_diags:          2

  arrowhead_resolution: 10000
  arrowhead_norm:       "SCALE"          # SCALE | KR
  arrowhead_java_mem:   "64g"

  hicexplorer:                 true      # third caller (hicFindTADs); shipped config = ON
  hicexplorer_resolution:      10000
  hicexplorer_min_depth_mult:  3         # min depth = 3 × res  (30 kb)
  hicexplorer_max_depth_mult:  10        # max depth = 10 × res (100 kb)
  hicexplorer_step_mult:       1         # step      = 1 × res  (10 kb)
  hicexplorer_fdr:             0.05      # --thresholdComparisons
  hicexplorer_delta:           0.01      # --delta

  agreement_slop_bp: 10000               # ± tolerance for "agreeing" boundaries
```

| Knob | Options |
|---|---|
| `insulation_primary_window` | **`200000`** *(200 kb at 10 kb bins — the standard)*. `100000` → sharper, more boundaries; `500000` → fewer, stronger. |
| `insulation_threshold` | **`Li`** *(cross-entropy threshold, robust on bimodal histograms)* / `Otsu` *(slightly more conservative)*. |
| `arrowhead_norm` | **`SCALE`** *(current Juicer recommendation)* / `KR` *(can return NaNs on sparse data)*. |
| `hicexplorer` | `true` adds hicFindTADs and switches the agreement metric to 3-caller mode; `false` → 2-caller mode (insulation + Arrowhead). |
| `agreement_slop_bp` | **`10000`** *(represents one 10 kb bin)*. `20000` for more lenient matching. |

### 3.5 Loop-caller settings

```yaml
loops:
  hiccups:     false                      # requires a CUDA GPU — off on BioQuant `single`
  chromosight: true                       # CPU, cooler-native — default here
  chromosight_resolution: 10000

  mustache:    false                      # opt-in second CPU caller → enables the loop-caller consensus
  mustache_resolution: 10000
  mustache_pthreshold: 0.1                # -pt (FDR)
  mustache_sthreshold: 0.88               # -st (sparsity)
  mustache_norm: "weight"                 # cooler's stored ICE weight column
  agreement_slop_bp: 10000                # both anchors must match within ± this interval

  # --- HiCCUPS / APA / HiCCUPSDiff (only used when hiccups: true) ---
  resolutions:    [5000, 10000]
  norm:           "SCALE"
  apa_resolution: 10000
  hiccups_java_mem: "96g"
  apa_java_mem:     "16g"
  diff_resolution:  10000
  diff_java_mem:    "32g"
  diff_pairs: []                          # HiCCUPSDiff contrasts, "A B" = loops of A relative to B,
                                          #   e.g. - "treated control"
```

| Knob | Options |
|---|---|
| `hiccups` | **`false`** on BioQuant `single` (no GPU). Switch to `true` only on a GPU cluster — Juicer compiles a CUDA kernel with `nvcc` at start-up. |
| `chromosight` | **`true`**. CPU-only, minutes per sample, often more sensitive than HiCCUPS at typical library depths. Only the resolution is exposed. |
| `mustache` | `true` runs Mustache as a second CPU caller; with ≥ 2 callers enabled the `loop_caller_agreement` rules run automatically (§7). |
| `diff_pairs` | HiCCUPSDiff contrasts, GPU only. |

### 3.6 Compartment settings

```yaml
compartments:
  resolution: 100000              # 100 kb standard; 250000 for shallow/sparse libraries
  saddle: true                    # compartment strength from a saddle plot
  n_quantile_bins: 38             # saddle quantile bins
  saddle_extreme_fraction: 0.2    # top/bottom 20 % of bins form the AA / BB / AB corners
```

| Knob | Options |
|---|---|
| `resolution` | **`100000`** *(Lieberman-Aiden standard)*. `250000` for shallow libraries. |
| `saddle` | `true` → `compartments/saddle/<sample>.saddle.txt` with `compartment_strength = mean(AA,BB) / mean(AB,BA)`. |
| `n_quantile_bins` / `saddle_extreme_fraction` | **`38`** / **`0.2`**; `0.1` gives a stricter strength metric. |

### 3.7 Optional modules

All optional modules are  OFF by default but can be switched ON in the config file; see §8 for details.

```yaml
runrule.trim:          false   # § 8.1 — fastp adapter-only trimming
loops.mustache:        false   # § 8.2 — second CPU loop caller + consensus (opt-in, not yet validated)
loops.hiccups:         false   # § 8.3 — HiCCUPS + APA + HiCCUPSDiff (GPU)
stripes.enabled:       false   # § 8.4 — Stripenn architectural stripes
differential.hicrep:   false   # § 8.5 — HiCRep pairwise similarity
differential.multihiccompare / .dchic: false   # § 8.5 — replicate-aware differential analysis
higlass.enabled:       false   # § 8.6 — HiGlass ingest manifest + launcher
```

---

## 4. Pipeline execution

Ideal start status: inside `tmux` on `login1`, Singularity loaded, `snakemake-9.19.0` env active, working
directory at the pipeline root, and a valid Kerberos ticket (see §1).

### 4.1 Dry-run

```bash
snakemake -n --configfile config/config.yaml
```

`-n` parses the sample sheet, validates every input path, builds the DAG and prints what would run (without
scheduling jobs or pulling containers). Read the metadata block at the top and the job-count summary at the
bottom. (`workflow/Snakefile` is found automatically; `--snakefile workflow/Snakefile` is optional.)

### 4.2 Launching the pipeline

```bash
snakemake --profile profiles/slurm \
          --configfile config/config.yaml \
          --rerun-incomplete --keep-going
```

- `--profile profiles/slurm` reads `profiles/slurm/config.yaml`, which sets the SLURM executor, apptainer-only
  software deployment, the bind mounts and the default resources (8 GB / 4 h / partition `single`).
- `--rerun-incomplete` redoes any rule whose outputs are flagged incomplete.
- `--keep-going` — when one sample's job fails, continue scheduling jobs for the other samples.

The first ever run additionally builds the bwa index (~1–2 h, one-time per genome) and pulls every container
image (~2.6 GB, one-time).

> **⚙ Adapt the SLURM profile to your cluster** — edit `profiles/slurm/config.yaml`: the `apptainer-args`
> bind mounts (every directory containing reads, genome, jar and results must be listed), `slurm-logdir`,
> `slurm_partition`, and the `slurm_extra="--exclude=…"` node list.

Per-rule resources (set inside the rules) for the heavy steps, per sample:

| Rule | threads | memory | wall-time |
|---|---|---|---|
| `align_parse_sort` (bwa + pairtools) | 16 | 64 GB | 48 h |
| `pairs_dedup` / `pairs_select` / `merge_lanes` | 4 | 16 GB | 8 h / 4 h / 4 h |
| `cooler_cload` → `cooler_zoomify` | 8 | 32 GB → 48 GB | 8 h → 12 h |
| `juicer_hic` (`juicer_tools pre`) | 8 | 80 GB | 24 h |
| `arrowhead` | 4 | 72 GB | 8 h |
| `insulation` / `hicfindtads` / `compartments_eigs` | 8 / 8 / 4 | 32 GB | 4 h / 8 h / 2 h |
| `chromosight_loops` / `mustache_loops` | 4 | 16 GB / 24 GB | 4 h / 6 h |
| `hiccups` (GPU only) | 8 | 110 GB | 12 h |

### 4.3 Monitoring progress

```bash
squeue -u "$USER"
tail -f .snakemake/log/$(ls -t .snakemake/log/*.snakemake.log | head -1)
ls results/logs/<stage>/          # per-rule tool logs; SLURM stdout/err under results/logs/slurm/rule_<name>/
```

### 4.4 Run locally without SLURM (small tests / debugging)

```bash
snakemake --cores 8 --software-deployment-method apptainer \
          --apptainer-args "--bind /path/to/reads --bind /path/to/genome --bind /path/to/results" \
          --configfile config/config.yaml
```

**Note**: alignment needs ~64 GB and `juicer_tools pre` ~80 GB of RAM per sample (only for tests).

### 4.5 Visualise the workflow

```bash
snakemake --rulegraph --configfile config/config.yaml | dot -Tpdf > rulegraph.pdf
snakemake --dag       --configfile config/config.yaml | dot -Tpdf > dag.pdf
```

---

## 5. Outputs & Interpretation

Under `<datapath>/`:

```
qc/multiqc/multiqc_report.html                 ← open this first (FastQC + pairtools dedup stats)
qc/fastqc/<sample>.<unit>_R{1,2}_fastqc.{zip,html}
qc/pairs_summary.tsv                           cross-sample: total / mapped / valid pairs / dup rate / cis / long-cis / trans
qc/scaling/<sample>.scaling.tsv                per-sample P(s) table (dist, n_valid, count_avg, balanced_avg)
qc/scaling/contact_decay.png                   overlay P(s) plot with a slope −1 reference

pairs/<sample>.valid.pairs.gz (+ .px2)         ★ the source of truth (UU/UR/RU, dedup, MAPQ ≥ 30, pairix-indexed)
pairs/<sample>.dedup.stats.txt                 pairtools dedup statistics (MultiQC-readable)

matrix/cool/<sample>.mcool                     ★ ICE-balanced pyramid 1 kb … 1 Mb (cooltools, chromosight, HiCExplorer, HiGlass)
matrix/hic/<sample>.hic                        ★ Juicer format with SCALE/KR (Arrowhead, HiCCUPS, Juicebox)

compartments/<sample>.cis.vecs.tsv             E1/E2/E3 per 100 kb bin (GC-phased)
compartments/<sample>.cis.lam.txt              eigenvalues per chromosome
compartments/<sample>.E1.bw                    E1 bigWig — load in IGV / HiGlass
compartments/saddle/<sample>.saddle.txt        compartment_strength scalar (+ .saddledump.npz matrix)
compartments/compartment_switching.tsv         per-bin A/B label per sample + `switches` flag

tads/insulation/<sample>.insulation.tsv        cooltools insulation table (all windows)
tads/insulation/<sample>.tads.bed              primary TAD set (segments between boundaries, 200 kb window)
tads/insulation/<sample>.boundaries.bed        boundary bins + strength
tads/insulation/boundary_dynamics.tsv          union of boundaries × sample presence (0/1) + n_samples
tads/arrowhead/<sample>.domains.bed            Juicer Arrowhead domains (BED3; raw 10000_blocks.bedpe kept in <sample>/)
tads/hicexplorer/<sample>.{domains,boundaries}.bed, <sample>.score.bedgraph   hicFindTADs (when tads.hicexplorer: true)
tads/agreement/<sample>.agreement.tsv          n per caller, pairwise Jaccards, consensus counts, frac-in-consensus
tads/agreement/<sample>.consensus.bed          BED5: boundaries called by ≥ 2 callers (name = callers, score = support)
tads/agreement/summary.tsv                     one row per sample

loops/chromosight/<sample>.chromosight.tsv     chromosight loops (+ .json, .pdf)
loops/mustache/<sample>.mustache.tsv           Mustache loops (when loops.mustache: true)
loops/agreement/<sample>.agreement.tsv         loop-caller Jaccards + consensus counts (≥ 2 callers enabled)
loops/agreement/<sample>.consensus.bedpe       ≥ 2-of-N consensus loops
loops/agreement/summary.tsv
loops/<sample>/merged_loops.bedpe, loops/<sample>/apa/, loops/diff/<A>_vs_<B>/    HiCCUPS / APA / HiCCUPSDiff (GPU only)

genome/analysis_view.bed, genome/gc_100000.tsv view + GC phasing track (built once)
logs/<stage>/*.log, logs/slurm/rule_<name>/    tool logs / SLURM job logs

# Optional modules (§8)
stripes/<sample>/result_{filtered,unfiltered}.tsv
differential/hicrep/<a>_vs_<b>.tsv, differential/hicrep/summary.tsv
differential/multihiccompare/<A>_vs_<B>.diff.tsv
differential/dchic/<A>_vs_<B>/dchic_{diff_compartments,AB_switches}.tsv
higlass/manifest.yaml, higlass/README.md, higlass/data/*.mcool
```

### 5.1 Library QC guide

| Metric | Source | Healthy characteristics (in mammalian Hi-C) |
|---|---|---|
| Duplicate rate | `qc/pairs_summary.tsv` (`dup_rate`), MultiQC | < 30–40 % |
| Valid pairs | `qc/pairs_summary.tsv` (`valid_pairs`) | more = better; ≥ 100 M for 10 kb loops |
| cis fraction (cis / (cis + trans)) | `qc/pairs_summary.tsv` (`cis_frac`) | ≥ 0.6–0.7; lower values = random ligation / low quality |
| long-cis fraction (≥ 20 kb) | `qc/pairs_summary.tsv` (`long_cis_frac_ge20kb`) | higher = better (≥ 0.4–0.5) |
| P(s) slope over 0.1–1 Mb | `qc/scaling/contact_decay.png` | ≈ −1 (log-log) |
| compartment strength | `compartments/saddle/<sample>.saddle.txt` | > 1; higher = cleaner A/B segregation |
| TAD-caller agreement | `tads/agreement/summary.tsv` (`frac_*_in_consensus`) | callers agree better if the libraries are clean and abundant |

**Note**: A library that fails on cis fraction, duplicate rate or P(s) most probably has a problem at the **library-level** (random ligation,
under-digestion, low complexity) -> cannot be fixed downstream. Sparse libraries typically give very few Arrowhead domains and chromosight loops while the insulation caller still reports boundaries.

---

## 6. TAD calling — three callers and the agreement metric

A "TAD" is not precisely defined - different callers legitimately disagree because they look for different
signal shapes (e.g., boundary-based vs domain-based; flat vs nested). Zufferey *et al.* (2018) show there is no gold
standard. The present pipeline therefore runs **three independent callers** and reports cross-caller agreement (consensus).

| Caller | Target | Output | Notes |
|---|---|---|---|
| **cooltools insulation** | local minima of the diamond-window insulation score (Crane *et al.* 2015), Li-threshold on boundary strength | `tads/insulation/*.tads.bed`, `*.boundaries.bed` | **primary** output, robust to resolution changes |
| **Juicer Arrowhead** | corner score of contact domains on the `.hic` (Rao *et al.* 2014) | `tads/arrowhead/*.domains.bed` | orthogonal, domain-based, handles nested TADs; needs abundant data |
| **HiCExplorer hicFindTADs** | insulation-style score swept over 30–100 kb diamonds with FDR-corrected boundary calls (Ramírez *et al.* 2018) | `tads/hicexplorer/*.domains.bed`, `*.boundaries.bed`, `*.score.bedgraph` | FDR-controlled alternative; `tads.hicexplorer: true` |

`rule tad_caller_agreement` (`workflow/scripts/tad_caller_agreement.py`):

- boundary points = midpoints of the boundary BEDs (insulation, hicFindTADs) and both endpoints of every
  Arrowhead domain;
- pairwise **slop-tolerant Jaccards** (`agreement_slop_bp`, default one 10 kb bin);
- single-linkage clustering of all boundary points across callers within the slop → each cluster's support =
  number of distinct callers;
- **consensus BED5** of boundaries supported by ≥ 2 callers; `summary.tsv` collects one row per sample
  (`n_<caller>`, `jaccard_*`, `n_consensus_3of3 / 2of3 / 1of3`, `frac_<caller>_in_consensus`).

With `hicexplorer: false`, the pipeline switches to 2-caller mode (`n_consensus_2of2 / 1of2`).
`boundary_dynamics.tsv` additionally clusters the insulation boundaries across samples (±10 kb) and marks in
which samples each boundary is present.

---

## 7. Loop calling — CPU (chromosight, Mustache) vs GPU (HiCCUPS)

| Caller | Input | GPU required? | Runtime / sample | Role in the present config |
|---|---|---|---|---|
| **chromosight** (Matthey-Doret *et al.* 2020) | `.mcool` | no | minutes | **default**|
| **Mustache** (Roayaei Ardakany *et al.* 2020) | `.mcool` | no | minutes–hour | 2nd CPU caller|
| **HiCCUPS** (Rao *et al.* 2014) | `.hic` | **yes** | hours | brings APA and HiCCUPSDiff when activated |

When ≥ 2 loop callers are enabled, `rule loop_caller_agreement` mirrors the TAD agreement in 2-D: anchors must
match within ± `loops.agreement_slop_bp` on **both** ends, pairwise Jaccards are reported and a ≥ 2-of-N
consensus BEDPE is written to `loops/agreement/`.

On a GPU cluster:  switch `hiccups: true`.

**Note**: HiCCUPS (Juicer's reference loop caller) **requires a CUDA GPU** — it compiles a PTX kernel via `nvcc` at
start-up and exits with `Cannot run program "nvcc"` / `GPU/CUDA Installation Not Detected` otherwise (error on the BioQuant `single` partition). For this reason, the present config uses loop callers that only use **CPU**.

---

## 8. Optional modules guide

- can be activated in the config file;
- OFF by default

| Module | Activation | Output | Status |
|---|---|---|---|
| **fastp trimming** | `runrule.trim: true` | `trim/<sample>.<unit>_R{1,2}.trimmed.fastq.gz` | not recommended for Hi-C |
| **Mustache + loop consensus** | `loops.mustache: true` | `loops/mustache/`, `loops/agreement/` | **currently untested** |
| **HiCCUPS / APA / HiCCUPSDiff** | `loops.hiccups: true` | `loops/<sample>/`, `loops/diff/` | needs a CUDA GPU; fails on BioQuant `single` partition|
| **Stripenn** | `stripes.enabled: true` | `stripes/<sample>/` | **currently untested** |
| **HiCRep** | `differential.hicrep: true` | `differential/hicrep/` | **currently untested** |
| **multiHiCcompare / differential E1** | `differential.multihiccompare` / `.dchic: true` + `condition` column | `differential/multihiccompare/`, `differential/dchic/` | **currently untested**; needs ≥ 2 conditions × ≥ 2 replicates |
| **HiGlass** | `higlass.enabled: true` | `higlass/` | **currently untested**; server needs Docker/Podman |

### 8.1 fastp adapter-only trimming

```yaml
runrule:
  trim: true
```
Runs `fastp` per lane with `--disable_quality_filtering --disable_length_filtering` (adapter removal only) and
feeds the trimmed FASTQs into `bwa mem`. Use is recommended only for documented adapter contamination (§3.2).

### 8.2 Mustache — 2nd CPU loop caller + loop-agreement consensus

```yaml
loops:
  mustache: true
  mustache_resolution: 10000
  mustache_pthreshold: 0.1        # -pt
  mustache_sthreshold: 0.88       # -st
  mustache_norm: "weight"         # cooler ICE weights (Mustache's default for .mcool)
  agreement_slop_bp: 10000
```

### 8.3 HiCCUPS, APA and HiCCUPSDiff

```yaml
loops:
  hiccups: true
  resolutions: [5000, 10000]
  norm: "SCALE"
  diff_pairs: ["A B", …]
```

### 8.4 Stripenn — architectural stripes (cohesin)

- detects asymmetric stripes that extend from a single anchor (~ snapshots of cohesin-mediated loop extrusion)

```yaml
stripes:
  enabled:   true
  resolution: 10000
  norm:       "weight"       # bins-table column with the balancing weights (cooler zoomify --balance writes "weight";
                          
  canny:      2.0
  max_pixel:  0.95
```

### 8.5 Replicate-aware differential analysis

```yaml
differential:
  hicrep:           true        # pairwise stratum-adjusted correlation (any pair of samples)
  hicrep_resolution: 50000
  hicrep_smoothing_h: 5
  hicrep_max_distance_bp: 5000000
  multihiccompare:  true        # differential bin pairs (needs ≥ 2 conditions × ≥ 2 reps)
  multihiccompare_resolution: 100000
  multihiccompare_fdr: 0.05
  dchic:            true        # differential compartments (same replicate gating)
  dchic_resolution: 100000
```

`samples.tsv` must carry `condition` (and `replicate`) columns for the last two. When no condition has ≥ 2
distinct samples, multiHiCcompare and the differential-compartment rule skip; 
HiCRep runs on every pair of samples regardless. 

The R modules install their packages on first use into a persistent `.R_libs/` at the
pipeline root (file-system lock serialises parallel installs) — the same pattern as the ChIP-Seq pipeline.


### 8.6 HiGlass — interactive web viewer

```yaml
higlass:
  enabled: true
  port:    8989
```

Then, on a host with Docker or Podman:

```bash
bash workflow/scripts/launch_higlass.sh \
    --data     <datapath>/higlass/data \
    --manifest <datapath>/higlass/manifest.yaml \
    --port 8989
```
Opens at `http://localhost:8989`.

---

## 9. Reproducibility — containers and provenance

Every rule declares a `container:` image and the SLURM profile only runs on **Apptainer**:

```yaml
software-deployment-method:
  - apptainer
apptainer-args: "--bind /net/bq-storage/ag-cherrmann/projects/36_NeuralCrest_Laugsch --bind /net/bq-storage/ag-cherrmann/projects/38_Pipelines"
rerun-triggers: "mtime"
slurm-logdir: ".../Hi-C/results/logs/slurm"
default-resources: [mem_mb=8000, runtime=240, slurm_partition=single, slurm_extra="--exclude=cln128,cln131,cln132"]
```

The bind mounts: project 36 (raw reads + Juicer jar) + project 38 (genome + pipeline + results) -> inside every container
Image registry (`workflow/rules/common.smk`):

| Pinned image | Used by |
|---|---|
| `docker://duplexa/4dn-hic:v43` | `align_parse_sort`, `merge_lanes`, `pairs_dedup`, `pairs_select`, `juicer_pairs` — bwa + pairtools + samtools + pairix + bgzip in one image (needed for the bwa→pairtools stream) |
| `docker://quay.io/biocontainers/bwa:0.7.18--h577a1d6_2` | `bwa_index` |
| `docker://quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3` | `gc_track`, `cooler_cload`, `cooler_zoomify`, `insulation`, `insulation_to_tads`, `boundary_dynamics`, `tad_caller_agreement`, `compartments_eigs`, `e1_bedgraph`, `compartment_saddle`, `compartment_switching`, `contact_decay(_plot)`, `pairs_summary`, `loop_caller_agreement`, `hicrep_pair` |
| `docker://eclipse-temurin:17` | `juicer_hic`, `arrowhead`, `hiccups`, `apa`, `hiccupsdiff` (Juicer Tools 2.20.00 jar) |
| `docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | `fastqc` |
| `docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5` | `fastp_trim` (optional) |
| `docker://ewels/multiqc:latest` | `multiqc` — the only **unpinned** tag |
| `docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:455--h2a80c09_1` | `e1_bigwig` |
| `docker://quay.io/biocontainers/chromosight:1.6.3--pyhdfd78af_0` | `chromosight_loops` |
| `docker://quay.io/biocontainers/hicexplorer:3.7.6--pyhdfd78af_0` | `hicfindtads` |
| `docker://quay.io/biocontainers/mustache-hic:1.3.3--pyhdfd78af_0` | `mustache_loops` |
| `docker://quay.io/biocontainers/stripenn:1.1.65.20--pyhdfd78af_0` | `stripenn_compute` |
| `docker://bioconductor/bioconductor_docker:RELEASE_3_19` | `multihiccompare_contrast`, `dchic_contrast` (+ `.R_libs/` install-on-first-use) |

Rules without a container (`analysis_view`, the three `*_summary` concatenations, `prep_higlass`) only use
Python/coreutils from the Snakemake environment.

The profile sets `rerun-triggers: mtime` — Snakemake re-runs jobs only when their **input files** change (by
mtime), not when a rule's code, parameters, or container reference changes. Routine actions such as editing comments or adding optional
modules therefore do not trigger a complete rerun of the pipeline.


## 10. References

### Workflow infrastructure
- **Snakemake.** Mölder F, *et al.* F1000Research 2021;10:33.
- **Apptainer / Singularity.** Kurtzer GM *et al.* PLoS ONE 2017;12:e0177459.

### Alignment + pair processing (4DN standard)
- **bwa-mem.** Li H. arXiv 2013;1303.3997.
- **pairtools.** Open2C *et al.* PLoS Comput Biol 2024;20:e1012164.
- **4DN Hi-C processing standard.** Reiff SB *et al.* Nat Commun 2022;13:2365.
- **pairix / 4dn-hic image.** 4D Nucleome DCIC, github.com/4dn-dcic/pairix, github.com/4dn-dcic/docker-4dn-hic.

### Matrix containers and balancing
- **cooler.** Abdennur N, Mirny LA. Bioinformatics 2020;36:311.
- **cooltools.** Open2C *et al.* PLoS Comput Biol 2024;20:e1012067.
- **ICE balancing.** Imakaev M *et al.* Nat Methods 2012;9:999.
- **KR balancing.** Knight P, Ruiz D. IMA J Numer Anal 2013;33:1029.
- **Juicer / Juicer Tools (SCALE, `pre`, Arrowhead, HiCCUPS, APA).** Durand NC *et al.* Cell Syst 2016;3:95; Rao SSP *et al.* Cell 2014;159:1665.

### A/B compartments
- **Compartments by eigendecomposition.** Lieberman-Aiden E *et al.* Science 2009;326:289.
- **Saddle plots / compartment strength.** Imakaev M *et al.* 2012 (above); cooltools documentation.

### TADs
- **Insulation score.** Crane E *et al.* Nature 2015;523:240.
- **Arrowhead.** Rao SSP *et al.* Cell 2014;159:1665.
- **HiCExplorer hicFindTADs.** Ramírez F *et al.* Nat Commun 2018;9:189.
- **TAD-caller benchmark.** Zufferey M *et al.* Genome Biol 2018;19:217.

### Loops and stripes
- **HiCCUPS / APA.** Rao SSP *et al.* Cell 2014;159:1665.
- **chromosight.** Matthey-Doret C *et al.* Nat Commun 2020;11:5795.
- **Mustache.** Roayaei Ardakany A *et al.* Genome Biol 2020;21:256.
- **Stripenn.** Yoon S *et al.* Nat Commun 2022;13:1602.

### Differential / reproducibility analysis
- **HiCRep.** Yang T *et al.* Genome Res 2017;27:1939.
- **multiHiCcompare.** Stansfield JC *et al.* Bioinformatics 2019;35:2916.
- **dcHiC.** Chakraborty A *et al.* Nat Commun 2022;13:6827.
- **limma.** Ritchie ME *et al.* NAR 2015;43:e47.

### QC / visualisation
- **FastQC.** Andrews S. babraham.ac.uk. **fastp.** Chen S *et al.* Bioinformatics 2018;34:i884.
- **MultiQC.** Ewels P *et al.* Bioinformatics 2016;32:3047.
- **HiGlass.** Kerpedjiev P *et al.* Genome Biol 2018;19:125.
- **Polymer physics of P(s).** Mirny LA. Chromosome Res 2011;19:37.
