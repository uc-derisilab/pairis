# PAIRIS: Prediction of Antibody-antigen Interactions at high Resolution In Silico

A Nextflow pipeline for high-throughput prediction of antibody-peptide complex structures using AlphaFold3 with optional Rosetta binding energy analysis, designed for SLURM HPC clusters.

## Features

- **MSA Reuse Optimization**: Pre-compute MSAs once and reuse across multiple complex predictions (4-6x speedup)
- **Sliding Window Analysis**: Generate overlapping peptide windows for epitope mapping
- **Template Integration**: Automatically extract and incorporate homology templates from AlphaFold3
- **Scalable Parallel Execution**: SLURM array jobs for efficient HPC scheduling
- **Modular Workflow**: Run individual stages independently or as a complete pipeline

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥24.10.0
- [AlphaFold3](https://github.com/google-deepmind/alphafold3) ≥3.0
- Python 3.12
- SLURM HPC cluster

**Python Dependencies:**

Main pipeline (MSA generation and structure prediction):
- No external packages required (uses Python standard library only)

Data collation:
- `pandas` ≥2.2.3

Rosetta analysis (optional):
- `pyrosetta` ≥2025.25 
- `biopython` ≥1.85 
- `polars` ≥1.32.3 

## Quick Start

### 1. Prepare Input Files

Create a project directory with your input sequences:

```bash
mkdir -p my_project/input/bcrs
cd my_project

# Add peptide sequences (one or more targets)
cat > input/peptides.fasta <<EOF
>Peptide_A
MKKIAVDKEGIPVALRR
>Peptide_B
EEVTGRGSQVEAFESREGGPWGGRVEAEESAGAEDSCGLDPAGSQTARA
EOF

# Add BCR sequences (one file per antibody, both chains)
cat > input/bcrs/mAb1.fasta <<EOF
>mAb1_heavy
QVQLVQSGAEVKKPGSSVKVSCKASGGTF...
>mAb1_light
DIQMTQSPSSLSASVGDRVTITCRASQSI...
EOF
```

**BCR Chain Requirements:**
- Each file must contain exactly two sequences (heavy and light chains)
- Headers must include `heavy`/`light` or `HC`/`LC` keywords (case-insensitive)

### 2. Configure Pipeline

Copy and edit the parameter template:

```bash
cp params.template.yml params.yml
vim params.yml
```

**Key parameters:**
```yaml
# Project settings
project_name: "my_project"
peptide_fasta: "input/peptides.fasta"
bcr_dir: "input/bcrs"

# Sliding window sizes (or null for full-length)
window_sizes: [15, 25]

# AlphaFold3 model diversity
num_seeds: 100

# Workflow stages
run_msa_generation: true
run_structure_prediction: true
run_rosetta_analysis: false

# Output directories
outdir: "results"
af3_output_dir: "/path/to/large_storage"  # For multi-GB AF3 outputs
```

### 3. Run Pipeline

```bash
nextflow run main.nf -params-file params.yml
```

**Resume interrupted runs:**
```bash
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

Predicts antibody-peptide complex structures using AlphaFold3:

- **Input:** Complex JSON files (auto-generated from Stage 1 or fresh MSAs)
- **MSA Mode:** Uses pre-computed MSAs if available (`--run_data_pipeline=false`)
- **Output:** Predicted structures (`.cif` files), confidence scores (`.json`)
- **Templates:** Automatically incorporated if found in MSA generation

### Stage 3: Rosetta Analysis (Optional)

Computes binding energies and interface metrics:

- **Input:** Predicted structures from Stage 2
- **Output:** CSV files with ∆G scores and interface residue analysis
- **Requirement:** Rosetta software installed

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

Override default resources in `params.yml`:

```yaml
# MSA generation (CPU-only)
msa_cpus: 2
msa_memory: '16.GB'
msa_time: '1.h'

# Structure folding (GPU-accelerated)
folding_cpus: 2
folding_memory: '16.GB'
folding_time: '40.m'
folding_gpu: 1

# Rosetta analysis (CPU-only)
rosetta_cpus: 4
rosetta_memory: '16.GB'
rosetta_time: '8.h'
```

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
