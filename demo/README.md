# PAIRIS demo

A small, self-contained example that predicts the structure of the **9E10
anti-c-myc antibody Fab bound to its c-myc epitope peptide** and (optionally)
reuses the precomputed MSAs. It exercises the full core pipeline — MSA
generation, complex input generation, and AlphaFold3 folding — on the smallest
realistic antibody–peptide system.

## Data provenance

All sequences are public and taken from **PDB [2OR9](https://www.rcsb.org/structure/2OR9)**
(murine monoclonal anti-c-myc antibody 9E10 in complex with its epitope peptide):

| File | Contents | PDB chain(s) |
|------|----------|--------------|
| `input/peptides.fasta` | c-myc epitope peptide `EQKLISEEDLN` | 2OR9 chain P |
| `input/bcrs/9E10.fasta` | 9E10 Fab heavy + light chains | 2OR9 chains H / L |

Because the deposited complex 2OR9 is available, you can compare the predicted
model against the experimental structure as an additional sanity check.

## Prerequisites

This demo **requires an NVIDIA GPU** plus a configured AlphaFold3 installation
(container `.sif`, model parameters, and genetic databases). See the
**System requirements** and **Installation** sections of the [top-level
README](../README.md) for how to obtain these. AlphaFold3 cannot run on a CPU-only
machine, so this demo cannot be run on a typical desktop without a GPU.

## Run the demo (single GPU server)

1. Edit `params.yml` and fill in the placeholder paths: `af3_output_dir`,
   `af3_sif`, `af3_model_dir`, `af3_db_dir`, and `local_venv`.
2. Create the output directory you chose for `af3_output_dir`:
   ```bash
   mkdir -p /path/to/demo_outputs
   ```
3. Launch from inside this directory:
   ```bash
   cd demo
   nextflow run ../main.nf -params-file params.yml -profile local
   ```

## Expected output

On success you will find (under the configured `af3_output_dir`):

```
af3_outputs/
├── msa_only/                       # MSA generation outputs (one per chain)
└── complexes/
    └── cmyc_epitope_9E10/
        └── cmyc_epitope_9e10/      # AF3 lowercases the job name
            ├── *_model.cif                  # predicted 3-chain complex structure
            ├── *_confidences.json           # per-residue/per-atom confidences
            └── *_summary_confidences.json   # overall iPTM / pLDDT scores
results/
└── msas/
    ├── msa_index.json              # sequence -> MSA path index
    └── msas/
        ├── cmyc_epitope/           # peptide MSA (shallow; short synthetic tag)
        ├── 9e10_heavy/             # heavy-chain MSA + templates
        └── 9e10_light/             # light-chain MSA + templates
```

The small per-run reports are written under `results/reports/`
(`report.html`, `timeline.html`, `trace.txt`). For reference, a validation run
predicted this complex with **iPTM ≈ pTM ≈ 0.81** (high-confidence).

## Expected run time

**~30 minutes end-to-end** for this configuration (validated on a single-GPU run).
MSA generation dominates: the three chains' database searches run in parallel and
take **~28 minutes** on CPU; AlphaFold3 folding on the GPU takes **~1.5 minutes**
for this ~450-residue system with `num_seeds: 1` / `num_diffusion_samples: 1`.
Enabling the optional Rosetta step adds under a minute. Actual time depends on
CPU/GPU hardware, database storage speed, and scheduler queueing; per-task timings
are recorded in `results/reports/trace.txt` and `report.html`.

## Alternative backend: ESMFold2

Two ready-to-edit configs run this same 9E10 / c-myc demo with ESMFold2 instead
of AlphaFold3 (see the [folding-backend docs](../README.md#folding-backend)):

- [`params.esmfold2.yml`](params.esmfold2.yml) — **ESMFold2-Fast**, MSA-free.
  Skips Stage 1 entirely and needs **no AlphaFold3 install or databases** — only
  an NVIDIA GPU, the `esmfold2` env, and the cached `biohub/ESMFold2-Fast` model:
  ```bash
  cd demo
  nextflow run ../main.nf -params-file params.esmfold2.yml -profile local
  ```
- [`params.esmfold2.full.yml`](params.esmfold2.full.yml) — **ESMFold2-full**,
  which uses MSAs. AlphaFold3's data pipeline is PAIRIS's only MSA source, so this
  still requires AlphaFold3 (it runs Stage 1) and folds with `biohub/ESMFold2`:
  ```bash
  cd demo
  nextflow run ../main.nf -params-file params.esmfold2.full.yml   # SLURM
  ```

Both write `*_model.cif` + `*_summary_confidences.json` in the same layout as the
AlphaFold3 backend, so the optional Rosetta + collation steps run unchanged.
Pre-download the model on a node with network access and point `esmfold2_hf_home`
at the cache.

## Quick check without a GPU (input generation only)

The pure-Python input-generation step uses only the standard library and runs in
under a second on any machine — useful to confirm the inputs parse correctly
before committing GPU time:

```bash
# from the repo root, using the demo data
python bin/generate_complex_inputs.py \
    --peptide-fasta demo/input/peptides.fasta \
    --antibody-fasta demo/input/bcrs/9E10.fasta \
    --antibody-name 9E10 \
    --num-seeds 1 \
    --output-dir /tmp/pairis_demo_inputs
```

This should report the heavy/light chain lengths and write a valid AlphaFold3
complex JSON (`cmyc_epitope_9E10.json`) — no GPU required.
