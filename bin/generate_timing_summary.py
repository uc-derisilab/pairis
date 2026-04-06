#!/usr/bin/env python3
"""Parse Nextflow trace.txt and produce a human-readable per-stage timing summary.

Requires trace.txt generated with `raw = true` in nextflow.config so that
duration/realtime columns are integer milliseconds.

Usage:
    generate_timing_summary.py --trace <trace.txt> --output <timing_summary.txt>
"""

import argparse
import csv
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path


STAGE_ORDER = [
    'GEN_MSA_PEPTIDES',
    'GEN_MSA_BCRS',
    'RUN_AF3_MSA',
    'RUN_ALPHAFAST_MSA',
    'EXTRACT_MSAS',
    'GENERATE_COMPLEX_INPUTS',
    'RUN_AF3_FOLDING',
    'GROUP_STRUCTURES_BY_BCR',
    'RUN_ROSETTA',
    'COLLATE_RESULTS',
]


def fmt_ms(ms):
    ms = int(ms)
    if ms < 1000:
        return f'{ms}ms'
    elif ms < 60_000:
        return f'{ms / 1000:.1f}s'
    elif ms < 3_600_000:
        return f'{ms / 60_000:.1f}min'
    else:
        return f'{ms / 3_600_000:.2f}h'


def stage_name(task_name):
    return re.sub(r'\s*\(.*\)$', '', task_name)


def parse_trace(trace_path):
    """Return dict: full_task_name -> realtime_ms (latest submit per name)."""
    tasks = {}  # full_name -> list of (submit, realtime_ms)
    with open(trace_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('status') != 'COMPLETED':
                continue
            name = row['name']
            try:
                realtime = int(row['realtime'])
            except (ValueError, KeyError):
                continue
            submit = row.get('submit', '')
            tasks.setdefault(name, []).append((submit, realtime))

    # Keep only the latest-submit run per task name
    final = {}
    for name, runs in tasks.items():
        final[name] = max(runs, key=lambda x: x[0])[1]
    return final


def build_summary(final_runs, pipeline_duration_str=None):
    stage_times = {}
    for full_name, rt in final_runs.items():
        s = stage_name(full_name)
        stage_times.setdefault(s, []).append(rt)

    lines = []
    lines.append('PAIRIS pipeline timing summary')
    lines.append(f'Generated:  {datetime.now().strftime("%Y-%m-%d %H:%M")}')
    if pipeline_duration_str:
        lines.append(f'Total wall clock: {pipeline_duration_str}')
    lines.append('')

    if not stage_times:
        lines.append('No COMPLETED tasks found (all tasks were CACHED or pipeline was skipped).')
        return '\n'.join(lines) + '\n'

    lines.append('=' * 72)
    lines.append('STAGE SUMMARY  (final successful run per task, CACHED excluded)')
    lines.append('=' * 72)
    lines.append(f'{"Stage":<30} {"N":>4}  {"Min":>8}  {"Median":>8}  {"Max":>8}  {"Total":>9}')
    lines.append('-' * 72)

    for s in STAGE_ORDER:
        if s not in stage_times:
            continue
        t = sorted(stage_times[s])
        n = len(t)
        med = statistics.median(t)
        lines.append(
            f'{s:<30} {n:>4}  {fmt_ms(min(t)):>8}  {fmt_ms(med):>8}'
            f'  {fmt_ms(max(t)):>8}  {fmt_ms(sum(t)):>9}'
        )
    # Any stages not in STAGE_ORDER (e.g. future additions)
    for s in sorted(stage_times):
        if s not in STAGE_ORDER:
            t = sorted(stage_times[s])
            n = len(t)
            med = statistics.median(t)
            lines.append(
                f'{s:<30} {n:>4}  {fmt_ms(min(t)):>8}  {fmt_ms(med):>8}'
                f'  {fmt_ms(max(t)):>8}  {fmt_ms(sum(t)):>9}'
            )

    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trace', required=True, help='Path to Nextflow trace.txt')
    parser.add_argument('--output', required=True, help='Path to write timing_summary.txt')
    parser.add_argument('--pipeline-duration', default=None,
                        help='Total pipeline wall-clock duration string (optional)')
    args = parser.parse_args()

    trace_path = Path(args.trace)
    if not trace_path.exists():
        print(f'ERROR: trace file not found: {trace_path}', file=sys.stderr)
        sys.exit(1)

    final_runs = parse_trace(trace_path)
    summary = build_summary(final_runs, args.pipeline_duration)

    Path(args.output).write_text(summary)
    print(summary, end='')


if __name__ == '__main__':
    main()
