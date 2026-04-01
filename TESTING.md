# PAIRIS Pipeline Testing Guide

This guide walks you through testing the PAIRIS pipeline step-by-step, focusing on validating the MSA fix and overall pipeline functionality.

## Conventions

Path placeholders used throughout this guide:
- `<project_dir>` — Your project working directory (e.g., where `params.yml` and `input/` live)
- `<af3_output_dir>` — The AF3 output directory set in `params.yml` (shared/group storage)
- `main.nf` — Run from the pipeline directory, or use the full path if running from elsewhere
- `bin/` — Scripts in the pipeline's `bin/` directory

## Overview

The pipeline has been updated to fix MSA file overwrites and optimize resource allocation. This testing guide will help you validate:

1. **MSA Generation & Extraction** - Ensures unique MSA directories are created
2. **Complex Input Generation** - Verifies MSA paths are inserted into JSONs
3. **Structure Prediction** - Validates AF3 uses precomputed MSAs
4. **Rosetta Analysis** - Confirms downstream analysis works

---

## Prerequisites

### 1. Project Directory Setup

Your project should have this structure:
```
<project_dir>/
├── params.yml              # Pipeline configuration
├── input/
│   ├── BCAS3_epitope.fasta
│   └── bcrs_test/
│       └── *.fasta
└── results/                # Will be created by pipeline
```

### 2. Verify params.yml

```bash
cd <project_dir>
cat params.yml
```

Should contain:
```yaml
project_name: "bcas3-mab"
peptide_fasta: "input/BCAS3_epitope.fasta"
bcr_dir: "input/bcrs_test"
window_sizes: [15]  # Start with just 15mer for testing
num_seeds: 5        # Reduced for quick testing
run_msa_generation: true
run_structure_prediction: true
run_rosetta_analysis: true
outdir: "results"
af3_output_dir: "<af3_output_dir>"

# SLURM partition settings (optional - defaults shown)
# partition: "preempted"     # Default for unlabeled processes
# msa_partition: "cpu"       # MSA generation (default: cpu)
# folding_partition: "gpu"   # Structure folding (default: gpu)
# rosetta_partition: "cpu"   # Rosetta analysis (default: cpu)
# extract_partition: "cpu"   # MSA extraction (default: cpu)
```

---

## Stage 1: MSA Generation & Extraction

### Goal
Validate that MSAs are extracted to unique directories without overwrites.

### Step 1.1: Clean Previous Results

```bash
cd <project_dir>

# Remove previous MSA results
rm -rf results/msas results/msa_index.json results/af3_inputs/msa

# Optional: Clean all previous results for fresh start
rm -rf work/ results/ .nextflow* <af3_output_dir>/af3_outputs/
```

### Step 1.2: Run MSA Generation Only

```bash
nextflow run main.nf \
  -params-file params.yml \
  --num_seeds 5 \
  --run_structure_prediction false \
  --run_rosetta_analysis false
```

**Note:** `window_sizes` is already set to `[15]` in params.yml. To override on command line, use comma-separated values:
```bash
# Single window size
--window_sizes 15

# Multiple window sizes
--window_sizes 15,25

# No sliding windows (full sequences)
--window_sizes null
```

**Expected behavior:**
- Jobs submit to `cpu` partition by default (can be overridden with `msa_partition` in params.yml)
- MSA generation doesn't need GPU, so `cpu` is recommended for faster queue times
- Pipeline runs MSA generation → MSA extraction → MSA index building
- Takes approximately 1-2 hours depending on queue

**To use a different partition (e.g., preempted):**
Add to your params.yml:
```yaml
msa_partition: "preempted"
extract_partition: "preempted"
```

Or override on command line:
```bash
nextflow run ... --msa_partition preempted --extract_partition preempted
```

### Step 1.3: Monitor Progress

```bash
# Watch the log
tail -f .nextflow.log

# Check SLURM queue
squeue -u $USER

# Check for errors
grep -i error .nextflow.log
```

### Step 1.4: Validate Results

After completion, run these validation checks:

```bash
# 1. Count unique MSA directories (should be 14: 12 peptide windows + heavy + light)
find results/msas/ -type d -mindepth 1 -maxdepth 1 | wc -l
# Expected output: 14

# 2. List all MSA directories to verify naming
find results/msas/ -type d -mindepth 1 -maxdepth 1 | sort
# Expected: BCAS3_epitope_window_000 through 011, mAb1_HA1C_heavy, mAb1_HA1C_light

# 3. Check MSA index exists and has 14 sequences
cat results/msa_index.json | jq 'keys | length'
# Expected output: 14

# 4. Verify each directory has required files
for dir in results/msas/*/; do
  echo "Checking $dir"
  ls -1 "$dir" | grep -E "(unpaired_msa.a3m|paired_msa.a3m|metadata.json)"
done
# Expected: Each directory should show all 3 files

# 5. Check for overwrite warnings (should see validation passed message)
grep -i "validation\|duplicate" .nextflow.log
# Expected: "Validation passed: 14 unique MSA directories created"
```

### ✅ Stage 1 Success Criteria

- [ ] 14 unique directories in `results/msas/`
- [ ] `msa_index.json` contains 14 entries
- [ ] No duplicate directory warnings in logs
- [ ] Each MSA directory has unpaired_msa.a3m, paired_msa.a3m, metadata.json

---

## Stage 2: Complex Input Generation

### Goal
Verify that complex JSONs are generated with MSA paths correctly inserted.

### Step 2.1: Generate Complex Inputs

```bash
nextflow run main.nf \
  -params-file params.yml \
  --num_seeds 5 \
  --run_msa_generation false \
  --run_structure_prediction true \
  --run_rosetta_analysis false \
  -resume
```

**Note:** This will generate complex JSONs but may also start AF3 folding. We'll validate JSONs before folding completes. The `window_sizes` setting from params.yml will be used.

### Step 2.2: Validate Complex JSONs

While the pipeline is running or after GENERATE_COMPLEX_INPUTS completes:

```bash
# 1. Check JSON version is 2 (required for MSA paths)
cat results/af3_inputs/complexes/*.json | jq -r '.version' | sort -u
# Expected output: 2

# 2. Check one complex JSON structure
cat results/af3_inputs/complexes/bcas3-mab_BCAS3_epitope_window_000_mAb1_HA1C.json | jq '.'
# Should show well-formed JSON with MSA paths

# 3. Verify MSA paths are present in all sequences
cat results/af3_inputs/complexes/*.json | jq '.sequences[].protein | select(has("unpairedMsaPath"))' | head -10
# Should show entries with unpairedMsaPath fields

# 4. Count sequences with MSA paths in each complex
for file in results/af3_inputs/complexes/*.json; do
  echo "File: $(basename $file)"
  jq '.sequences | map(select(.protein.unpairedMsaPath != null)) | length' "$file"
done
# Expected: 3 for each file (peptide + heavy + light chains)

# 5. Validate all MSA paths exist on filesystem
jq -r '.sequences[].protein.unpairedMsaPath // empty' results/af3_inputs/complexes/*.json | while read path; do
  if [ ! -f "$path" ]; then
    echo "MISSING: $path"
  fi
done
# Expected: No output (all paths exist)
```

### ✅ Stage 2 Success Criteria

- [ ] All complex JSONs have `version: 2`
- [ ] All 3 sequences in each complex have `unpairedMsaPath` fields
- [ ] All MSA paths point to existing files
- [ ] JSON structure is valid

---

## Stage 3: Structure Prediction

### Goal
Validate that AF3 uses precomputed MSAs and completes successfully.

### Step 3.1: Run Structure Prediction

If you stopped the pipeline after Stage 2, restart with:

```bash
nextflow run main.nf \
  -params-file params.yml \
  --num_seeds 5 \
  --run_msa_generation false \
  --run_structure_prediction true \
  --run_rosetta_analysis false \
  -resume
```

**Expected behavior:**
- Jobs submit to `gpu` partition by default (can be overridden with `folding_partition` in params.yml)
- AF3 folding requires GPU, so `gpu` partition is necessary
- Pipeline uses precomputed MSAs from Stage 1
- Uses `window_sizes` from params.yml
- Faster than without precomputed MSAs (~40 min per structure vs ~2-3h)

**Note:** If your cluster has a different GPU partition name, override in params.yml:
```yaml
folding_partition: "gpu-a100"  # or whatever your GPU partition is called
```

### Step 3.2: Monitor Folding Progress

```bash
# Check running jobs
squeue -u $USER | grep gpu

# Monitor specific job logs
tail -f work/<job_hash>/.command.log

# Check overall progress
grep "RUN_AF3_FOLDING" .nextflow.log | tail -20
```

### Step 3.3: Validate Structure Outputs

```bash
# 1. Check AF3 output directories exist
ls -d <af3_output_dir>/af3_outputs/complexes/*/
# Should show multiple directories

# 2. Count structure files generated
find <af3_output_dir>/af3_outputs/complexes/ -name "*_model_*.cif" | wc -l
# Expected: At least 1 per complex (depends on num_seeds)

# 3. Check a structure file exists
ls -lh <af3_output_dir>/af3_outputs/complexes/bcas3-mab_BCAS3_epitope_window_000_mAb1_HA1C/
# Should see *_model_*.cif files

# 4. Verify cleanup happened (seed folders removed)
find <af3_output_dir>/af3_outputs/complexes/ -name "seed_*" -type d | wc -l
# Expected: 0 (all cleaned up)

# 5. Check execution time from Nextflow logs
grep "RUN_AF3_FOLDING" .nextflow.log | grep -E "COMPLETED" | head -5
# Note the duration - should be faster than runs without MSAs
```

### ✅ Stage 3 Success Criteria

- [ ] Structure prediction completes without errors
- [ ] CIF files generated for all complexes
- [ ] Seed folders are cleaned up (saves disk space)
- [ ] Runtime is reduced compared to baseline

---

## Stage 4: Rosetta Analysis (Optional)

### Goal
Validate that Rosetta analysis works with generated structures.

### Step 4.1: Run Rosetta Analysis

```bash
nextflow run main.nf \
  -params-file params.yml \
  --run_msa_generation false \
  --run_structure_prediction false \
  --run_rosetta_analysis true \
  --rosetta_time '1.h' \
  -resume
```

**Note:** Using `--rosetta_time '1.h'` for quick testing (default is 8h for full datasets)

### Step 4.2: Validate Rosetta Results

```bash
# 1. Check Rosetta CSV files generated
ls -lh results/rosetta/*.csv | wc -l

# 2. View a sample CSV
head -20 results/rosetta/BCAS3_epitope_window_000_rosetta.csv

# 3. Check combined results file
head -20 results/reports/combined_results.csv

# 4. View summary statistics
cat results/reports/summary_stats.txt
```

### ✅ Stage 4 Success Criteria

- [ ] Rosetta CSV files generated for all structures
- [ ] Combined results file exists
- [ ] Summary statistics computed

---

## Full Pipeline Test

### Goal
Run complete pipeline end-to-end to validate all stages work together.

### Step 5.1: Clean Slate

```bash
cd <project_dir>
rm -rf work/ results/ .nextflow* <af3_output_dir>/af3_outputs/
```

### Step 5.2: Run Full Pipeline

```bash
nextflow run main.nf \
  -params-file params.yml \
  --num_seeds 5
```

**Expected behavior:**
- MSA generation runs on `cpu` partition
- Structure folding runs on `gpu` partition
- Rosetta runs on `cpu` partition
- Uses `window_sizes: [15]` from params.yml
- All stages complete in sequence

### Step 5.3: Final Validation

Run all validation checks from Stages 1-4 above to ensure:
- MSAs extracted correctly
- Complex JSONs have MSA paths
- Structures predicted successfully
- Rosetta analysis completes

---

## Expected Performance

Based on the MSA fix:

| Metric | Before Fix | After Fix |
|--------|------------|-----------|
| MSA files created | 2 (12 overwrites) | 14 (unique) |
| MSA reuse | None | All complexes |
| Structure prediction time | ~2-3h per complex | ~40min per complex |
| Total pipeline runtime | ~24h | ~12h (~50% reduction) |

---

## Troubleshooting

### Issue: MSAs still overwriting

**Symptoms:** `find results/msas/ -type d | wc -l` returns less than 14

**Solution:**
```bash
# Check the extraction script
grep "json_file.stem.replace" bin/extract_msas_from_precomputed.py
# Should show: seq_id = json_file.stem.replace('_data', '')

# Check logs for errors
grep -A5 "extract_msas_from_precomputed" .nextflow.log
```

### Issue: Complex JSONs missing MSA paths

**Symptoms:** `jq '.sequences[0].protein | has("unpairedMsaPath")' file.json` returns false

**Solution:**
```bash
# Check if MSA index was found
grep "msa_index" .nextflow.log

# Verify MSA index exists
cat results/msa_index.json | jq 'keys'

# Check GENERATE_COMPLEX_INPUTS received the index
grep "GENERATE_COMPLEX_INPUTS" .nextflow.log | grep "msa"
```

### Issue: Jobs stuck in queue

**Symptoms:** `squeue -u $USER` shows jobs pending for long time

**Solution:**
```bash
# Check partition availability
sinfo -p cpu
sinfo -p gpu

# Try different partition
nextflow run ... --partition preempted

# Reduce resource requirements
nextflow run ... --msa_cpus 2 --msa_memory '8.GB'
```

### Issue: JSON version is still 1

**Symptoms:** `jq '.version' file.json` returns 1 instead of 2

**Solution:**
```bash
# Check utils.py has version 2
grep "version.*2" bin/utils.py

# If not, update the file or contact pipeline maintainer
```

### Issue: Pipeline fails with "permission denied"

**Solution:**
```bash
# Make sure scripts are executable
chmod +x bin/*.py

# Check output directory permissions
ls -ld results/
ls -ld <af3_output_dir>/
```

---

## Quick Reference Commands

### Resume failed pipeline
```bash
nextflow run main.nf \
  -params-file params.yml \
  -resume
```

### Override resources for faster testing
```bash
nextflow run main.nf \
  -params-file params.yml \
  --num_seeds 1 \
  --window_sizes 15 \
  --folding_time '20.m'
```

**Note:** Use comma-separated values for multiple window sizes: `--window_sizes 15,25`

### Run only specific stage
```bash
# MSA only
nextflow run ... --run_structure_prediction false --run_rosetta_analysis false

# Structure only (requires existing MSAs)
nextflow run ... --run_msa_generation false --run_rosetta_analysis false

# Rosetta only (requires existing structures)
nextflow run ... --run_msa_generation false --run_structure_prediction false
```

### Check pipeline status
```bash
# View execution report (after completion)
firefox results/reports/report.html

# View timeline
firefox results/reports/timeline.html

# Check detailed trace
less results/reports/trace.txt
```

---

## Getting Help

### Log Files Locations

- **Nextflow log:** `.nextflow.log` (in project directory)
- **Process logs:** `work/<process_hash>/.command.log`
- **SLURM logs:** Check `work/<process_hash>/.command.err`
- **Execution reports:** `results/reports/`

### Useful Debug Commands

```bash
# Show all failed processes
grep "ERROR" .nextflow.log

# Show specific process details
nextflow log -f name,status,duration,workdir

# Get work directory for specific process
nextflow log -f name,workdir | grep RUN_AF3_FOLDING
```

### Contact

For issues with:
- **Pipeline bugs:** Check existing code or contact maintainer
- **SLURM/cluster issues:** Contact HPC support
- **AlphaFold3 errors:** Check AF3 documentation
- **Nextflow errors:** Check Nextflow documentation

---

## Summary Checklist

After completing all testing stages, you should have:

- [ ] 14 unique MSA directories in `results/msas/`
- [ ] MSA index with 14 entries in `results/msa_index.json`
- [ ] Complex JSONs with version 2 and MSA paths
- [ ] Structure CIF files in AF3 output directory
- [ ] Rosetta CSV files in `results/rosetta/`
- [ ] Execution reports in `results/reports/`
- [ ] No errors in `.nextflow.log`
- [ ] Faster runtime compared to baseline

**If all checkboxes are checked, the pipeline is working correctly!** ✅
