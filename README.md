# PAIRIS: Prediction of Antibody-antigen Interactions at high Resolution In Silico

A Nextflow pipeline for high-throughput prediction of antibody–peptide complex structures,
with optional Rosetta binding-energy analysis. Structure prediction can use **AlphaFold3**
(default) or **ESMFold2** — a fast, optionally MSA-free alternative — selectable via the
`folding_backend` parameter. PAIRIS runs on SLURM HPC clusters (default) or on a single GPU
server (`-profile local`).

## Features

- **Selectable Folding Backend**: AlphaFold3 (default) or ESMFold2 (a fast, optionally MSA-free alternative) via `folding_backend`
- **MSA Reuse Optimization**: Pre-compute MSAs once and reuse across multiple complex predictions (4-6x speedup)
- **Sliding Window Analysis**: Generate overlapping peptide windows for epitope mapping
- **Template Integration**: Automatically extract and incorporate homology templates from AlphaFold3
- **Scalable Parallel Execution**: SLURM array jobs for efficient HPC scheduling
- **Modular Workflow**: Run individual stages independently or as a complete pipeline

## System requirements

### Operating system

- **Linux, x86_64 only** (AlphaFold3/Apptainer are Linux-only; no macOS/Windows).
- **Tested on** Rocky Linux 8.10, kernel `4.18.0-553.el8`.

### Software dependencies

PAIRIS is pure Nextflow + Python (no compilation). The external tools it orchestrates must be
installed separately:

| Component | Required | Tested version |
|-----------|----------|----------------|
| [Nextflow](https://www.nextflow.io/) | ≥24.10.0 | 24.10.5 |
| Java (Nextflow runtime) | 17+ | — |
| [Apptainer](https://apptainer.org/) / Singularity (for `-profile local`) | any recent | 1.4.5 |
| [AlphaFold3](https://github.com/google-deepmind/alphafold3) | ≥3.0 (default backend; also the only MSA source) | 3.0.1 (module `alphafold/3.0.1-23-g792e61e`) |
| ESMFold2 (`esm` Python package; alternative backend, `folding_backend: esmfold2`) | optional | see `env/requirements-esmfold2.txt` |
| Python | 3.12 | 3.12.7 |
| SLURM (only for the default cluster profile) | — | — |

> **AlphaFold3 model parameters and genetic databases are not redistributed here.** Request
> the parameters from Google DeepMind and download the databases separately — see
> [Installation](#installation) and [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

### Python packages

Input/MSA generation uses **only the Python standard library**. Two small environments (conda
or `uv`) cover the rest; name them anything and point `conda_env` / `rosetta_conda_env` at them:

- **Default** — [`env/requirements.txt`](env/requirements.txt), for result collation. Python
  3.12.7, `pandas` 2.2.3 (the only hard dependency; `polars` 0.20.31 and `biopython` 1.84 are
  pinned to match the tested env).
- **Rosetta** (optional, when `run_rosetta_analysis: true`) —
  [`env/requirements-rosetta.txt`](env/requirements-rosetta.txt). Python 3.12.11, `polars`
  1.32.3, `biopython` 1.85, `pyrosetta` 2025.25 (not on PyPI; needs a Rosetta license — see
  that file and THIRD_PARTY_LICENSES.md).
- **ESMFold2** (optional, when `folding_backend: esmfold2`) —
  [`env/requirements-esmfold2.txt`](env/requirements-esmfold2.txt), used by
  `bin/run_esmfold2.py`; point `esmfold2_conda_env` at it (default `esmfold2`). Python 3.12,
  `torch`, `pydantic`, and a forked `esm` that transitively installs a forked `transformers`
  providing `transformers.models.esmfold2` (plain PyPI `transformers` will not work).
  **Requires an NVIDIA GPU**; the model is downloaded from HuggingFace on first run (set
  `esmfold2_hf_home` to a cache dir and pre-download on a node with network access).

### Hardware (non-standard)

- **NVIDIA GPU required** for inference with either backend — typically a data-center GPU
  (A100/H100-class, ~40–80 GB); small complexes like the demo run on less. *(Tested on an
  NVIDIA A100 80GB.)*
- **~hundreds of GB of disk** for the AF3 genetic databases. The ESMFold2 backend with
  `esmfold2_model: fast` is MSA-free and needs no AF3 databases at all (just the HuggingFace
  model cache); `esmfold2_model: full` still needs AF3-generated MSAs.
- A normal desktop **cannot** run AlphaFold3 or ESMFold2 inference, though the
  input-generation steps run on any CPU (see the [demo](#demo)).

## Installation

No compilation required. Times below assume a normal desktop/login node with good network.

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd pairis
   ```
2. **Install Nextflow** (≥24.10.0) and Java 17+.
3. **Install Apptainer/Singularity** for `-profile local` (often preinstalled on HPC).
4. **Set up AlphaFold3** (one-time, several hours): obtain the container (`.sif`), request model
   parameters from Google DeepMind, and download the genetic databases (the multi-hundred-GB
   download dominates). See the [AlphaFold3 guide](https://github.com/google-deepmind/alphafold3).
5. **Create the Python environment(s)**:
   ```bash
   # default environment (any name; referenced via conda_env, default "pairis")
   conda create -n pairis python=3.12
   conda activate pairis
   pip install -r env/requirements.txt

   # optional: Rosetta analysis environment
   conda create -n rosetta python=3.12
   conda activate rosetta
   pip install -r env/requirements-rosetta.txt
   # then install pyrosetta per env/requirements-rosetta.txt
   ```
6. **(Optional) Create the ESMFold2 environment** — only needed for
   `folding_backend: esmfold2`. Use a separate env so its forked `esm`/`transformers` do not
   clash with the default env (`mamba` works too; the name is referenced via
   `esmfold2_conda_env`):
   ```bash
   conda create -n esmfold2 python=3.12
   conda activate esmfold2
   pip install -r env/requirements-esmfold2.txt
   ```
   The model is downloaded from HuggingFace Hub on first run. HPC compute nodes often have no
   internet, so set `esmfold2_hf_home` to a cache directory and pre-download the model on a
   node with network access (`biohub/ESMFold2-Fast` for the `fast` variant,
   `biohub/ESMFold2` for `full`).

**Typical install time:** PAIRIS + Python environment, **~10 minutes**. One-time AlphaFold3
setup (container + parameters + databases), **several hours** (dominated by the database
download) — a standard AF3 prerequisite, not specific to PAIRIS.

## Demo

[`demo/`](demo/) holds a small self-contained example — the **9E10 anti-c-myc Fab** bound to
its **c-myc epitope** (sequences from PDB 2OR9) — exercising the full core pipeline (MSA
generation, complex inputs, AF3 folding).

```bash
cd demo
# edit params.yml to fill in af3_output_dir, af3_sif, af3_model_dir, af3_db_dir
nextflow run ../main.nf -params-file params.yml -profile local
```

**Output:** a 3-chain complex structure (`*_model.cif`), confidence scores
(`*_summary_confidences.json`, iPTM/pLDDT ≈ 0.81 on a validation run), and MSAs with an
`msa_index.json`. **Run time:** ~**30 minutes** (a single-GPU run; dominated by the ~28 min
CPU MSA database search, with ~1.5 min GPU folding). The demo needs a GPU; full instructions
and a GPU-free sanity check are in [`demo/README.md`](demo/README.md).

## Instructions for use (your own data)

### 1. Prepare input files

```bash
mkdir -p my_project/input/bcrs
cd my_project

# Peptide / epitope sequences (one or more targets)
cat > input/peptides.fasta <<EOF
>Peptide_A
MKKIAVDKEGIPVALRR
>Peptide_B
EEVTGRGSQVEAFESREGGPWGGRVEAEESAGAEDSCGLDPAGSQTARA
EOF

# Antibody sequences: one FASTA per antibody, with both chains
cat > input/bcrs/mAb1.fasta <<EOF
>mAb1_heavy
QVQLVQSGAEVKKPGSSVKVSCKASGGTF...
>mAb1_light
DIQMTQSPSSLSASVGDRVTITCRASQSI...
EOF
```

**BCR chain requirements:**
- Each file must contain exactly two sequences (heavy and light chains).
- Headers must include `heavy`/`light` or `HC`/`LC` keywords (case-insensitive).

### 2. Configure the pipeline

```bash
cp params.template.yml params.yml
$EDITOR params.yml
```

Key parameters:
```yaml
project_name: "my_project"
peptide_fasta: "input/peptides.fasta"
bcr_dir: "input/bcrs"

window_sizes: [15, 25]   # or null for full-length peptides
num_seeds: 100           # AlphaFold3 model diversity

run_msa_generation: true
run_structure_prediction: true
run_rosetta_analysis: false

outdir: "results"
af3_output_dir: "/path/to/large_storage"   # multi-GB AF3 outputs
```

### 3. Run the pipeline

```bash
# On a SLURM cluster (default)
nextflow run main.nf -params-file params.yml

# On a single GPU server (fill in af3_sif / af3_model_dir / af3_db_dir first)
nextflow run main.nf -params-file params.yml -profile local

# Fold with ESMFold2-Fast (MSA-free; no AlphaFold3 needed)
nextflow run main.nf -params-file params.yml \
    --folding_backend esmfold2 --esmfold2_model fast --run_msa_generation false

# Resume an interrupted run
nextflow run main.nf -params-file params.yml -resume
```

## Pipeline Stages

### Stage 1: MSA Generation (Optional but Recommended)

Pre-computes multiple sequence alignments for all input sequences:

- **Input:** Peptide and BCR FASTA files
- **Output:** Extracted MSA files (`.a3m` format) + index
- **Benefit:** Reuse MSAs across all complex predictions (4-6x faster)

**Skip if:** You're running a single prediction and don't need optimization.

### Stage 2: Structure Prediction

Predicts antibody-peptide complex structures. The backend is selectable via `folding_backend`:

- **AlphaFold3** (default):
  - **Input:** Complex JSON files (auto-generated from Stage 1 or fresh MSAs)
  - **MSA Mode:** Uses pre-computed MSAs if available (`--run_data_pipeline=false`)
  - **Output:** Predicted structures (`.cif` files), confidence scores (`.json`)
  - **Templates:** Automatically incorporated if found in MSA generation
- **ESMFold2** (`folding_backend: esmfold2`): drop-in alternative that writes to the same
  output directories, so Rosetta and collation run unchanged. With `esmfold2_model: fast`
  it is MSA-free, so **Stage 1 can be skipped** and AlphaFold3 is not needed at all. With
  `esmfold2_model: full` it still consumes AF3-generated MSAs, so Stage 1 (or an existing
  `msa_index.json`) is required. For each complex it folds `num_seeds` seeds and keeps the
  single best by peptide-interface confidence.

### Stage 3: Rosetta Analysis (Optional)

Computes binding energies and interface metrics:

- **Input:** Predicted structures from Stage 2
- **Output:** CSV files with ∆G scores and interface residue analysis
- **Requirement:** Rosetta environment (see [System requirements](#python-packages))

## Output Structure

```
my_project/
├── results/
│   ├── af3_inputs/
│   │   ├── msa/              # MSA-only JSON inputs
│   │   └── complexes/        # BCR-peptide complex JSONs
│   ├── rosetta/              # Binding energy CSVs (if enabled)
│   └── reports/
│       ├── timeline.html     # Execution timeline
│       ├── report.html       # Resource usage
│       └── trace.txt         # Detailed logs
│
└── af3_outputs/              # Large outputs (stored separately)
    ├── msa_only/             # AlphaFold3 MSA generation outputs
    │   └── <seq_id>/
    │       └── <seed>/
    │           └── *_data.json
    ├── complexes/            # AlphaFold3 structure predictions
    │   └── <complex_id>/
    │       └── <seed>/
    │           ├── *_model.cif
    │           ├── *_confidences.json
    │           └── *_summary_confidences.json
    └── results/
        └── msas/
            ├── msa_index.json      # MSA path index
            └── msas/
                └── <seq_id>/
                    ├── unpaired_msa.a3m
                    ├── paired_msa.a3m
                    ├── metadata.json
                    └── templates/  # Homology templates (if found)
```

## Configuration

### Folding backend

`folding_backend` selects the structure-prediction engine:

```yaml
folding_backend: "alphafold3"   # default, unchanged behavior
# folding_backend: "esmfold2"   # fast, optionally MSA-free alternative
```

When `folding_backend: esmfold2`, these additional parameters apply:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `esmfold2_model` | `fast` | `fast` = ESMFold2-Fast, MSA-free (skips Stage 1 entirely; **no AlphaFold3 needed**). `full` = uses MSAs; since AF3's data pipeline is PAIRIS's only MSA source, `full` **still requires AlphaFold3** plus `run_msa_generation: true` or an existing `msa_index.json`. |
| `esmfold2_conda_env` | `esmfold2` | conda/mamba env activated for the ESMFold2 folding step (see [`env/requirements-esmfold2.txt`](env/requirements-esmfold2.txt)). |
| `esmfold2_hf_home` | `null` | HuggingFace cache dir (`HF_HOME`) for the downloaded model; set it and pre-download on a node with network access. |
| `esmfold2_num_loops` | `3` | ESMFold2 `num_loops`. |
| `esmfold2_num_sampling_steps` | `50` | ESMFold2 `num_sampling_steps`. |

ESMFold2 outputs land in the same `af3_inputs/complexes` and `af3_outputs/complexes`
directories as AlphaFold3, so the Rosetta and collation stages work unchanged.

### Sliding Window Analysis

Control peptide fragmentation for epitope mapping:

```yaml
# Both 15-mer and 25-mer overlapping windows
window_sizes: [15, 25]

# Single window size
window_sizes: [15]

# No sliding window (use full peptide sequence)
window_sizes: null
```

### Resource Tuning

Override default resources in `params.yml` (see [RESOURCES.md](RESOURCES.md) for the full table):

```yaml
# MSA generation (CPU-only)
msa_cpus: 2
msa_memory: '16.GB'
msa_time: '1.h'

# Structure folding (GPU-accelerated)
folding_cpus: 1
folding_memory: '16.GB'
folding_time: '40.m'
folding_gpu: 1

# Rosetta analysis (CPU-only)
rosetta_cpus: 4
rosetta_memory: '16.GB'
rosetta_time: '8.h'
```

> **ESMFold2 seeds run serially.** Unlike AlphaFold3's single parallel pass, the ESMFold2
> backend folds `num_seeds` seeds one at a time per complex (keeping the best by
> peptide-interface confidence). The ESMFold2 folding step reuses the `folding_*` resources,
> so with a large `num_seeds` (e.g. the default `100`) raise `folding_time` accordingly — or
> lower `num_seeds` — to avoid wall-clock timeouts.

### Independent Stage Execution

**MSA generation only:**
```bash
nextflow run main.nf -params-file params.yml \
    --run_structure_prediction false \
    --run_rosetta_analysis false
```

**Structure prediction using existing MSAs:**
```bash
nextflow run main.nf -params-file params.yml \
    --run_msa_generation false \
    --run_rosetta_analysis false
```

**Rosetta analysis on existing structures:**
```bash
nextflow run main.nf -params-file params.yml \
    --run_msa_generation false \
    --run_structure_prediction false
```

## Performance Optimizations

### MSA Reuse Strategy

The pipeline builds an MSA index during Stage 1, allowing Stage 2 to:
1. Skip redundant MSA computations for the same sequences
2. Load pre-computed MSAs directly via file paths
3. Incorporate homology templates automatically

**Expected speedup:** 4-6x faster structure prediction vs. fresh MSA generation

### SLURM Array Jobs

Large task batches are submitted as array jobs rather than individual jobs:
- Cleaner queue management
- Better scheduler efficiency
- Configured via `array = 10000` in `nextflow.config`

## Troubleshooting

### BCR Chain Detection Failure

**Error:** `Could not identify heavy and light chains in <file>`

**Solution:** Ensure chain headers contain required keywords:
```fasta
>mAb1_heavy  ✓
>mAb1_HC     ✓
>mAb1_chain1 ✗ (missing keyword)
```

### Missing MSA Paths in Complex JSONs

**Symptom:** Structure prediction starts MSA generation despite `run_msa_generation=false`

**Solution:** Verify MSA index exists and is correctly loaded:
```bash
# Check index file
cat <af3_output_dir>/results/msas/msa_index.json

# Ensure msa_index_file is set correctly in main.nf
```

## Further documentation

- [`demo/README.md`](demo/README.md) — runnable demo with expected output and run time
- [RESOURCES.md](RESOURCES.md) — per-process resource allocation and tuning
- [TESTING.md](TESTING.md) — stage-by-stage validation guide
- [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) — AlphaFold3 and Rosetta license terms
- [LICENSE](LICENSE) — MIT license for the PAIRIS source code
