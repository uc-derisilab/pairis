# Checking a PAIRIS Run

How to confirm each stage worked, and what to do when it didn't. Run the checks from your
project folder (where `params.yml` lives). `<af3_output_dir>` and `<outdir>` are the values in
your `params.yml` (`outdir` defaults to `results`).

**Tip:** for a first test of a new setup, use one antibody, one short peptide, `num_seeds: 1`
and `window_sizes: [15]`, or just run the [demo](demo/README.md).

## How many of everything to expect

| Thing | Expected count |
|-------|----------------|
| Peptide sequences to fold | per peptide: 1 if `window_sizes: null`, else `length − k + 1` windows for each window size `k` |
| MSA folders (stage 1) | peptide sequences + 2 per antibody (heavy, light) |
| Complexes (stage 2) | peptide sequences × antibodies |
| Collated CSVs (stage 3) | 1 per antibody |

Example: one 26-residue peptide with `window_sizes: [15]` gives 12 windows; with one antibody
that's 14 MSA folders and 12 complexes.

## Stage 1: MSAs

```bash
# MSA folders (one per sequence)
ls <af3_output_dir>/results/msas/msas/ | wc -l

# Index used by stage 2: should list the same number of entries
jq 'keys | length' <af3_output_dir>/results/msas/msa_index.json

# Every folder should have these three files
for d in <af3_output_dir>/results/msas/msas/*/; do
  for f in unpaired_msa.a3m paired_msa.a3m metadata.json; do
    [ -f "$d$f" ] || echo "MISSING: $d$f"
  done
done
```

**Good:** counts match the table above, and the loop prints nothing.

## Stage 2: Folding

**Inputs.** Each complex JSON should point every chain at an MSA file that exists:

```bash
# Every chain in every complex should have an MSA path (prints nothing if OK)
for f in <outdir>/af3_inputs/complexes/*.json; do
  n=$(jq '[.sequences[].protein | select(.unpairedMsaPath != null)] | length' "$f")
  [ "$n" = 3 ] || echo "$f: only $n of 3 chains have MSAs"
done

# Every MSA path should exist (prints nothing if OK)
jq -r '.sequences[].protein.unpairedMsaPath // empty' <outdir>/af3_inputs/complexes/*.json \
  | sort -u | while read p; do [ -f "$p" ] || echo "MISSING: $p"; done
```

These JSON checks apply to AlphaFold3. ESMFold2-Fast inputs have no MSAs.

**Outputs.** One structure per complex:

```bash
find <af3_output_dir>/af3_outputs/complexes -name "*_model.cif" | wc -l
```

**Good:** the count equals the number of complexes. If some are missing, re-run with `-resume`;
completed complexes are skipped.

## Stage 3: Rosetta and collation

```bash
# Per-peptide Rosetta results
find <af3_output_dir>/results/rosetta -name "*_rosetta.csv" | wc -l

# Final results: one CSV per antibody, plus a summary
ls <af3_output_dir>/results/collated/
cat <af3_output_dir>/results/collated/summary.txt
head <af3_output_dir>/results/collated/<antibody>.csv
```

**Good:** one collated CSV per antibody, with `iptm` and `binding_energy` filled in for each
peptide or window. Empty `binding_energy` values mean Rosetta failed for that peptide; see below.
For what the numbers mean, see [Outputs in the README](README.md#outputs).

## When something fails

**Find the failed job** (same in both modes):
```bash
grep -i error .nextflow.log | tail                                 # what Nextflow saw
nextflow log last -f name,status,exit,workdir | grep -v COMPLETED  # failed tasks + work folders
```
Then in that work folder, read `.command.err` and `.command.log` (the job's own output) and
`.command.sh` (the exact command it ran).

| Symptom | Likely cause and fix |
|---------|----------------------|
| Exit code 137/140, or "time limit" in the log | Ran out of memory or time. Raise `*_memory` / `*_time` ([RESOURCES.md](RESOURCES.md)). Retries alone won't help. |
| Fewer MSA folders than expected | Check the RUN_AF3_MSA logs; a sequence may have failed every retry. Re-run with `-resume`. |
| Complex JSONs missing MSA paths | `msa_index.json` not found: check `af3_output_dir` points at the run that made the MSAs. |
| Antibody missing from results | Its FASTA headers must contain `heavy`/`light` or `HC`/`LC`. |
| SLURM: jobs pending for a long time | Partition busy or wrong name. Check with `sinfo`, then set `*_partition` in `params.yml`. |
| Single server: `no free GPU slot` | A killed run left GPU locks behind. With no PAIRIS run active, `rmdir /tmp/pairis_gpu_*.lock`. Keep `local_max_parallel_folding` ≤ `local_num_gpus`. |
| Single server: runs slowly or the machine is overloaded | Lower `local_max_parallel_cpu`; MSA searches are CPU- and disk-heavy. Check GPU use with `nvidia-smi`. |
| Want a clean re-run of one stage | Delete only that stage's outputs (see the README's [Outputs](README.md#outputs)) and re-run with `-resume`. Don't delete `<af3_output_dir>` wholesale; it holds the expensive MSAs and structures. |

**Run reports** in `<outdir>/reports/`: `report.html` (time and memory per job), `timeline.html`
(when each job ran) and `trace.txt` (one line per job, with status and exit code).
