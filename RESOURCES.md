# PAIRIS Compute Resources

Default CPUs, memory, time limits and SLURM partitions for each step, and how to change them.
All defaults live in `nextflow.config`; override any `*_cpus`/`*_memory`/`*_time`/`*_partition`
value in `params.yml` or on the command line (e.g. `--folding_time 2.h`).

## Defaults

| Stage | Step | CPUs | Memory | Time | Partition | GPU | Retries |
|-------|------|------|--------|------|-----------|-----|---------|
| 1. MSAs | GENERATE_MSA_INPUTS | 1 | 4 GB | 5m | `msa_partition` (cpu) | - | - |
| | RUN_AF3_MSA | `msa_cpus` 2 | `msa_memory` 16 GB | `msa_time` 1h | `msa_partition` (cpu) | - | 3 |
| | EXTRACT_MSAS | `extract_cpus` 2 | `extract_memory` 16 GB | `extract_time` 60m | `extract_partition` (cpu) | - | - |
| 2. Folding | GENERATE_COMPLEX_INPUTS / GENERATE_ESMFOLD2_INPUTS | 1 | 4 GB | 5m | `msa_partition` (cpu) | - | - |
| | RUN_AF3_FOLDING / RUN_ESMFOLD2_FOLDING | `folding_cpus` 1 | `folding_memory` 16 GB | `folding_time` 40m | `folding_partition` (gpu) | `folding_gpu` 1 | 10 |
| 3. Rosetta | GROUP_STRUCTURES_BY_BCR | 1 | 8 GB | 2h | `grouping_partition` (cpu) | - | - |
| | RUN_ROSETTA | `rosetta_cpus` 4 | `rosetta_memory` 16 GB | `rosetta_time` 8h | `rosetta_partition` (cpu) | - | 10 |
| | COLLATE_RESULTS | 1 | 4 GB | 45m | `partition` (preempted) | - | - |

Names in `code` are the parameters you can override. Retries resubmit a failed job with the
**same** resources, so they help with transient failures (node preemption, storage hiccups) but
not with jobs that run out of memory or time. For those, raise the limit.

## Common adjustments

```yaml
# Folding timing out (large complexes, or ESMFold2 with many seeds, which run one after another)
folding_time: '2.h'

# Rosetta running out of memory
rosetta_memory: '32.GB'

# Your cluster uses different partition names
msa_partition: 'cpu'
folding_partition: 'gpu'
rosetta_partition: 'cpu'
partition: 'preempted'   # anything not listed above
```

## How jobs are submitted

- **Array jobs:** the big steps (MSAs, folding, Rosetta, input generation) are submitted as
  SLURM array jobs: one array with one task per sequence, complex or Rosetta group,
  instead of thousands of separate jobs.
- **Queue limits:** Nextflow keeps at most 500 jobs queued and submits at most 10 per second
  (`executor` block in `nextflow.config`).
- **Environments:** most steps run `module load mamba && mamba activate <conda_env>` (default
  `pairis`). RUN_ROSETTA uses `rosetta_conda_env` and ESMFold2 folding uses `esmfold2_conda_env`.
  AlphaFold3 steps load the module set by `af3_module` (default `alphafold/3.0.1-23-g792e61e`).
  With `-profile local` this differs: see below.

## Running on a single server (`-profile local`)

There is no SLURM, so partitions and arrays are ignored. Instead, these limit how many jobs run
at once:

| Parameter | Default | Limits |
|-----------|---------|--------|
| `local_max_parallel_folding` | 4 | concurrent GPU folding jobs and AF3 MSA jobs |
| `local_max_parallel_cpu` | 8 | concurrent Rosetta and input-generation jobs |
| `max_cpus` / `max_memory` | 32 / 128 GB | total the pipeline may use on the machine |

Environments in single-server mode: the general steps activate the venv at `local_venv`, and
AlphaFold3 runs from the `af3_sif` container. Rosetta and ESMFold2 activate the same conda envs
as on SLURM, but through `conda shell.bash hook` rather than `module load`, so `conda` must be
on your `PATH`.
