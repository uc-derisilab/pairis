# PAIRIS Pipeline Default Resources

This document describes the default computational resources allocated to each process in the PAIRIS pipeline.

## Process Resource Allocation

| Stage | Process | CPUs | Memory | Time | Partition | GPU | Array Size |
|-------|---------|------|--------|------|-----------|-----|------------|
| **MSA Generation** | GENERATE_MSA_INPUTS | 1 | 4 GB | 30m | preempted | - | - |
| | RUN_AF3_MSA | 2 | 16 GB | 1h | cpu | - | 10000 |
| | EXTRACT_MSAS | 2 | 16 GB | 45m | cpu | - | - |
| **Structure Prediction** | GENERATE_COMPLEX_INPUTS | 1 | 4 GB | 30m | preempted | - | - |
| | RUN_AF3_FOLDING | 2 | 16 GB | 40m | gpu | 1 | 10000 |
| **Rosetta Analysis** | GROUP_STRUCTURES_BY_BCR | 1 | 4 GB | 30m | preempted | - | - |
| | RUN_ROSETTA | 4 | 16 GB | 8h | cpu | - | 10000 |
| | COLLATE_RESULTS | 1 | 4 GB | 30m | preempted | - | - |

## Global Resource Limits

- **Max Memory**: 128 GB
- **Max CPUs**: 32
- **Max Time**: 24 hours
- **Default Partition**: preempted
- **SLURM Queue Size**: 500
- **Submit Rate Limit**: 10 submissions per second

## Process Details

### MSA Generation Stage

**GENERATE_MSA_INPUTS**
- Generates AlphaFold3 MSA input JSON files
- Lightweight process that creates input files only
- No special labels

**RUN_AF3_MSA** (label: `af3_msa`)
- Runs AlphaFold3 MSA generation only (no inference)
- Uses array submission for parallel processing (up to 10,000 jobs)
- Retry strategy: 3 attempts on failure

**EXTRACT_MSAS** (label: `extract`)
- Extracts and organizes MSA files from AF3 output
- Creates MSA index for complex generation

### Structure Prediction Stage

**GENERATE_COMPLEX_INPUTS**
- Generates AlphaFold3 complex input JSON files
- Creates peptide-antibody complex inputs with MSA references
- Lightweight process

**RUN_AF3_FOLDING** (label: `af3_folding`)
- Runs AlphaFold3 structure prediction on GPU
- Uses precomputed MSAs (run_data_pipeline=false)
- Requires 1 GPU via `--gres=gpu:1`
- Array submission for massive parallelization (up to 10,000 jobs)
- Retry strategy: 3 attempts on failure
- Automatically cleans up large intermediate files

### Rosetta Analysis Stage

**GROUP_STRUCTURES_BY_BCR**
- Groups AF3 structures by (kmer_size, BCR, peptide) for analysis
- Lightweight process that creates grouping metadata

**RUN_ROSETTA** (label: `rosetta`)
- Performs Rosetta InterfaceAnalyzer energy calculations
- Most resource-intensive CPU process (8 hour time limit)
- Uses the `rosetta_conda_env` environment instead of the default `conda_env`
- Array submission for parallel processing (up to 10,000 jobs)
- Retry strategy: 2 attempts on failure

**COLLATE_RESULTS**
- Combines AF3 metrics (iptm, peptide sequences) with Rosetta energies
- Lightweight process that creates final summary CSVs

## Customizing Resources

Resources can be overridden in `params.yml` or via command line:

### MSA Generation
```yaml
msa_cpus: 2
msa_memory: '16.GB'
msa_time: '1.h'
msa_partition: 'cpu'
```

### Structure Folding
```yaml
folding_cpus: 2
folding_memory: '16.GB'
folding_time: '40.m'
folding_gpu: 1
folding_partition: 'gpu'
```

### Rosetta Analysis
```yaml
rosetta_cpus: 4
rosetta_memory: '16.GB'
rosetta_time: '8.h'
rosetta_partition: 'cpu'
```

### Extraction
```yaml
extract_cpus: 2
extract_memory: '16.GB'
extract_time: '45.m'
extract_partition: 'cpu'
```

## Environment Configuration

- **Default Environment**: All processes run `mamba activate <conda_env>` (default `pairis`)
- **Rosetta Environment**: RUN_ROSETTA runs `mamba activate <rosetta_conda_env>` (default `rosetta`)
- **AlphaFold3**: RUN_AF3_MSA and RUN_AF3_FOLDING load `module load alphafold/3.0.1-23-g792e61e`

## Array Job Limits

Processes with `array: 10000` can submit up to 10,000 parallel jobs:
- RUN_AF3_MSA: One job per peptide sequence/window
- RUN_AF3_FOLDING: One job per peptide-antibody complex
- RUN_ROSETTA: One job per (kmer_size, BCR, peptide) combination

This allows for massive scalability when analyzing hundreds of BCRs or peptide sequences.

## Notes

1. **GPU Requirement**: Only RUN_AF3_FOLDING requires GPU resources
2. **Retry Logic**: AF3 processes retry up to 3 times, Rosetta retries up to 2 times
3. **Partition Strategy**:
   - `preempted`: Lightweight processes (generation, grouping, collation)
   - `cpu`: CPU-intensive processes (MSA, extraction, Rosetta)
   - `gpu`: GPU-required processes (folding)
4. **Time Allocation**: Rosetta analysis has the longest time limit (8 hours) due to energy calculations
5. **Resource Configuration**: See `nextflow.config` for full configuration details
