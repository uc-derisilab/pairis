# PAIRIS: Prediction of Antibody-antigen Interactions at high Resolution In Silico

PAIRIS predicts how antibodies bind peptides. You give it antibody sequences and peptide
sequences; it folds every antibody–peptide pair and tells you which peptides (or which
windows of a longer peptide) each antibody most likely binds.

It runs in up to three stages, each of which can be switched on or off:

| Stage | What it does | Tool |
|-------|--------------|------|
| 1. MSAs | Builds sequence alignments once per sequence, reused by every pair | AlphaFold3 (CPU) |
| 2. Folding | Predicts the 3D structure of each antibody–peptide complex | AlphaFold3 (default) or ESMFold2 (GPU) |
| 3. Rosetta | Scores binding energy of each predicted complex, then collates all results | PyRosetta (CPU) |

PAIRIS is a [Nextflow](https://www.nextflow.io/) pipeline. It runs in one of two modes:

| Mode | Use when | How it runs |
|------|----------|-------------|
| **Single server** (`-profile local`) | You have one Linux machine with NVIDIA GPU(s) | Jobs run directly on that machine; AlphaFold3 runs from an Apptainer container |
| **SLURM cluster** (default) | You have an HPC cluster with SLURM | Jobs are submitted to the queue; AlphaFold3 is loaded as a cluster module |

Everything else (inputs, settings, outputs) is the same in both modes.

## Quick start

```bash
# 1. Make a project folder with your inputs (see "Inputs" below)
mkdir -p my_project/input/bcrs && cd my_project

# 2. Copy and edit the settings file
cp /path/to/pairis/params.template.yml params.yml

# 3. Run inside tmux/screen: closing the terminal stops the pipeline
tmux new -s pairis

#    single server:
nextflow run /path/to/pairis/main.nf -params-file params.yml -profile local
#    SLURM cluster:
nextflow run /path/to/pairis/main.nf -params-file params.yml

# If it stops partway, re-run the same command with -resume to pick up where it left off
nextflow run /path/to/pairis/main.nf -params-file params.yml -profile local -resume
```

On SLURM you can instead submit Nextflow itself as a long-running job (e.g. an `sbatch` script
that runs the command above), so it keeps going after you log out.

The answer you want ends up in `<af3_output_dir>/results/collated/` (see [Outputs](#outputs)).

## Instructions for use

### Inputs

**Peptides:** one FASTA file, one or more sequences:
```
>Peptide_A
MKKIAVDKEGIPVALRR
>Peptide_B
EEVTGRGSQVEAFESREGGPWGGRVEAEESAGAEDSCGLDPAGSQTARA
```

**Antibodies:** one FASTA file per antibody in a folder (e.g. `input/bcrs/mAb1.fasta`),
each with exactly two sequences, heavy and light. Headers must contain `heavy`/`light` or
`HC`/`LC` (any case):
```
>mAb1_heavy
QVQLVQSGAEVKKPGSSVKVSCKASGGTF...
>mAb1_light
DIQMTQSPSSLSASVGDRVTITCRASQSI...
```

### Settings (`params.yml`)

`params.template.yml` lists every option with comments. The ones you will usually change:

```yaml
project_name: "my_project"
peptide_fasta: "input/peptides.fasta"
bcr_dir: "input/bcrs"
outdir: "results"                         # small files: inputs, run reports
af3_output_dir: "/path/to/large_storage"  # big files: structures, MSAs, final results

window_sizes: [15, 25]   # split peptides into overlapping 15- and 25-mers; null = whole peptide
num_seeds: 100           # predictions per complex; more = slower but more robust

run_msa_generation: true        # stage 1
run_structure_prediction: true  # stage 2
run_rosetta_analysis: true      # stage 3

folding_backend: "alphafold3"   # or "esmfold2" (see below)
```

**Single-server settings** (needed for `-profile local`):

```yaml
af3_sif: "/path/to/alphafold.sif"            # AlphaFold3 Apptainer container
af3_model_dir: "/path/to/af3/params"         # AlphaFold3 model weights
af3_db_dir: "/path/to/af3/public_databases"  # AlphaFold3 genetic databases
local_venv: "/path/to/pairis-venv"           # Python env for the general steps
local_num_gpus: 1                            # GPUs on this machine
local_max_parallel_folding: 1                # folding jobs at once; keep ≤ local_num_gpus
local_max_parallel_cpu: 8                    # CPU jobs at once (Rosetta, input generation)
```

The `af3_*` paths are only checked when AlphaFold3 will run, so ESMFold2-Fast runs can skip them.

**SLURM settings:** usually only partition names, if your cluster's differ from `cpu`/`gpu`/
`preempted` (see [RESOURCES.md](RESOURCES.md)), and `af3_module` if your AlphaFold3 module has
a different name.

**Sliding windows** are how PAIRIS maps epitopes: a long peptide is cut into overlapping
windows, each folded with the antibody, so you can see which region scores best.

### Running only some stages

Each stage reads the previous stage's outputs from `af3_output_dir`, so you can run them
separately, for example to re-score existing structures without re-folding:

```bash
# Stage 1 only
nextflow run main.nf -params-file params.yml --run_structure_prediction false --run_rosetta_analysis false

# Stage 2 only, reusing existing MSAs
nextflow run main.nf -params-file params.yml --run_msa_generation false --run_rosetta_analysis false

# Stage 3 only, on existing structures
nextflow run main.nf -params-file params.yml --run_msa_generation false --run_structure_prediction false
```

### Using ESMFold2 instead of AlphaFold3

ESMFold2 is a faster alternative for stage 2. Set `folding_backend: "esmfold2"` and pick a model:

| `esmfold2_model` | Needs MSAs? | Needs AlphaFold3? |
|------------------|-------------|-------------------|
| `fast` (default) | No: set `run_msa_generation: false` | No |
| `full` | Yes: from stage 1 or an earlier run | Yes, for stage 1 |

Other ESMFold2 settings: `esmfold2_conda_env` (default `esmfold2`), `esmfold2_hf_home` (model
cache folder; compute nodes usually have no internet, so download the model there first),
`esmfold2_num_loops` (3), `esmfold2_num_sampling_steps` (50).

ESMFold2 runs seeds one after another, so with `num_seeds: 100` raise `folding_time` or lower
`num_seeds` to avoid timeouts. It writes to the same folders as AlphaFold3, so stage 3 works
unchanged.

### Compute resources

Defaults (CPUs, memory, time, SLURM partition per stage) are in `nextflow.config` and can
be overridden in `params.yml`, e.g. `folding_time: '2.h'`. See [RESOURCES.md](RESOURCES.md).

## Outputs

```
<af3_output_dir>/
├── results/
│   ├── collated/                 ← START HERE
│   │   ├── <antibody>.csv        # one row per peptide/window: iptm, binding_energy
│   │   └── summary.txt
│   ├── rosetta/<k>/<antibody>/   # full Rosetta metrics per peptide
│   └── msas/                     # extracted MSAs + msa_index.json
└── af3_outputs/
    ├── msa_only/                 # raw stage 1 output
    └── complexes/<complex>/      # predicted structures (*_model.cif) + confidence JSONs

<outdir>/                          # default: results/ in your project folder
├── af3_inputs/                   # generated input JSONs
└── reports/                      # timeline.html, report.html, trace.txt (run diagnostics)
```

**Reading the collated CSVs** (columns: `peptide_name`, `window`, `peptide_sequence`,
`kmer_size`, `iptm`, `binding_energy`):
- `iptm`: AlphaFold3's confidence in the peptide's interface with the antibody (the peptide
  chain's ipTM), 0–1. **Higher is
  better**; roughly, above 0.8 is confident and below 0.6 is unreliable.
- `binding_energy`: Rosetta interface energy (Rosetta energy units). **More negative is
  stronger binding.** Compare it between peptides for the same antibody, not as an absolute
  affinity.

A peptide or window that scores well on both is the best epitope candidate.

## Demo

[`demo/`](demo/) is a small worked example: the 9E10 anti-c-myc antibody with its c-myc
epitope (PDB 2OR9). It takes about 30 minutes on one GPU and gives ipTM ≈ 0.81.

```bash
cd demo
# edit params.yml: fill in af3_output_dir, af3_sif, af3_model_dir, af3_db_dir, local_venv
nextflow run ../main.nf -params-file params.yml -profile local
```

See [`demo/README.md`](demo/README.md) for expected output and a check that runs without a GPU.

## System requirements

- **Linux x86_64** (tested on Rocky Linux 8.10)
- **NVIDIA GPU** for folding (tested on an A100 80GB; small complexes need less)
- **Disk:** hundreds of GB for the AlphaFold3 databases, plus space for outputs
  (`af3_output_dir`); large runs produce many GB

| Software | Version (tested) | Needed for |
|----------|------------------|------------|
| Nextflow + Java 17+ | ≥24.10.0 (24.10.5) | everything |
| Apptainer/Singularity | any recent (1.4.5) | single server: runs AlphaFold3 |
| uv (or any venv tool) | any | single server: Python env for the general steps |
| SLURM + environment modules | — | cluster mode |
| conda/mamba | — | Rosetta and ESMFold2 steps, in both modes: on a single server `conda` must be on your `PATH`; on SLURM it comes from `module load mamba` |
| AlphaFold3 | ≥3.0 (3.0.1) | stages 1–2 (not needed for ESMFold2 `fast`): a `.sif` container on a single server, a module on SLURM |
| Python | 3.12 | collation, Rosetta, ESMFold2 |
| PyRosetta | 2025.25 | stage 3 (needs a Rosetta license) |

AlphaFold3 model parameters and databases are **not** included: request the parameters from
Google DeepMind. See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

## Installation

About 10 minutes for PAIRIS. AlphaFold3 setup is a one-time, several-hour job, mostly the
database download.

1. Clone this repository.
2. Install Nextflow (≥24.10.0) and Java 17+.
3. Set up AlphaFold3 (skip this for ESMFold2-Fast only):
   - **Single server:** follow the [AlphaFold3 guide](https://github.com/google-deepmind/alphafold3)
     to get the `.sif` container, model parameters (request from Google DeepMind) and
     databases, then point `af3_sif`, `af3_model_dir`, `af3_db_dir` at them.
   - **SLURM:** the cluster's AlphaFold3 module supplies weights and databases. Set
     `af3_module` to its name (default `alphafold/3.0.1-23-g792e61e`).
4. Create the Python environments.

   General steps, **single server**: a venv, pointed to by `local_venv`:
   ```bash
   uv venv /path/to/pairis-venv && source /path/to/pairis-venv/bin/activate
   uv pip install -r env/requirements.txt
   ```
   General steps, **SLURM**: a conda env, named by `conda_env` (default `pairis`):
   ```bash
   conda create -n pairis python=3.12 && conda activate pairis
   pip install -r env/requirements.txt
   ```
   Rosetta and ESMFold2, **both modes**: conda envs, named by `rosetta_conda_env` and
   `esmfold2_conda_env`:
   ```bash
   # stage 3 (Rosetta): pyrosetta install steps are in the file
   conda create -n rosetta python=3.12 && conda activate rosetta
   pip install -r env/requirements-rosetta.txt

   # ESMFold2 backend only; keep it separate, it uses forked esm/transformers
   conda create -n esmfold2 python=3.12 && conda activate esmfold2
   pip install -r env/requirements-esmfold2.txt
   ```

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `Could not identify heavy and light chains in <file>` | Antibody FASTA headers must contain `heavy`/`light` or `HC`/`LC`; `>mAb1_chain1` won't work. |
| Stage 2 regenerates MSAs even with `run_msa_generation: false` | Check `<af3_output_dir>/results/msas/msa_index.json` exists; `af3_output_dir` must point at the run that made the MSAs. |
| Jobs hit their time limit | Raise that stage's `*_time` in `params.yml` (see [RESOURCES.md](RESOURCES.md)). |
| Single server: `no free GPU slot` | A killed run left GPU locks behind. With no PAIRIS run active, `rmdir /tmp/pairis_gpu_*.lock`. Also check `local_max_parallel_folding` ≤ `local_num_gpus`. |
| Single server: `conda: command not found` or `Could not find conda environment` in a Rosetta/ESMFold2 job | Make sure `conda` is on your `PATH` in the shell that launches Nextflow, or set `rosetta_conda_env` / `esmfold2_conda_env` to the env's full path (from `conda env list`). |
| Single server: `-profile local requires --af3_sif` | Fill in `af3_sif`, `af3_model_dir`, `af3_db_dir` in `params.yml`. |
| A job failed and you want details | Find the task's `hash` (e.g. `3a/1b2c4d`) in `<outdir>/reports/trace.txt`, then read `.command.log` and `.command.err` in `work/3a/1b2c4d*/`. |
| Run stopped partway | Re-run the same command with `-resume`. |

## More documentation

- [`params.template.yml`](params.template.yml): every setting, with comments
- [`demo/README.md`](demo/README.md): the runnable demo
- [RESOURCES.md](RESOURCES.md): resource defaults and tuning
- [TESTING.md](TESTING.md): checking each stage's output
- [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md): AlphaFold3 and Rosetta licenses
- [LICENSE](LICENSE): MIT license for PAIRIS
