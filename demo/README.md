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
        └── <seed>/
            ├── *_model.cif                  # predicted 3-chain complex structure
            ├── *_confidences.json           # per-residue/per-atom confidences
            └── *_summary_confidences.json   # overall iPTM / pLDDT scores
results/
└── msas/
    ├── msa_index.json              # sequence -> MSA path index
    └── msas/
        ├── cmyc_epitope/           # peptide MSA (shallow; short synthetic tag)
        ├── 9E10_heavy/             # heavy-chain MSA + templates
        └── 9E10_light/             # light-chain MSA + templates
```

The small per-run reports are written under `results/reports/`
(`report.html`, `timeline.html`, `trace.txt`).

## Expected run time

**Estimate: roughly 10–20 minutes** on a single modern data-center GPU
(e.g. NVIDIA A100/H100), most of it spent searching the genetic databases during
MSA generation; the folding step itself is short for this ~450-residue system
with `num_seeds: 1` and `num_diffusion_samples: 1`. Actual time depends on GPU
model, database storage speed, and scheduler overhead. After running, the
measured wall-clock is recorded in `results/reports/report.html`.

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
