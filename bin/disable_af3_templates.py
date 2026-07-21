#!/usr/bin/env python3
"""Rewrite an AlphaFold3 input JSON with template search disabled.

Sets ``"templates": []`` on every protein chain, which makes AF3 skip the
template hmmsearch entirely (run_template_search=False) while STILL running the
MSA/jackhmmer search. ``unpairedMsa``/``pairedMsa`` are left untouched (unset),
which is the code path AF3 requires for "search MSA, no templates".

Used by the RUN_AF3_MSA retry fallback: a task whose template search hangs past
its wall clock is retried on a JSON produced by this script.

Usage: disable_af3_templates.py <input.json> <output.json>
"""
import json
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: disable_af3_templates.py <input.json> <output.json>")
    src, dst = sys.argv[1], sys.argv[2]

    with open(src) as fh:
        data = json.load(fh)

    for entry in data.get("sequences", []):
        if "protein" in entry:
            # Empty list => AF3 sets run_template_search=False (skips template
            # hmmsearch). MSA fields are intentionally left unset so jackhmmer
            # still runs.
            entry["protein"]["templates"] = []

    with open(dst, "w") as fh:
        json.dump(data, fh, indent=4)


if __name__ == "__main__":
    main()
