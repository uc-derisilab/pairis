#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files for MSA generation (--run_inference=false).
Creates JSON files with protein sequences and a single model seed.
"""

import argparse
import json
import os


def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences as dictionary."""
    sequences = {}
    current_id = None

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line

    return sequences


def generate_sliding_windows(sequence, window_size):
    """Generate sliding windows of specified size across sequence."""
    if len(sequence) < window_size:
        return [sequence]

    windows = []
    for i in range(len(sequence) - window_size + 1):
        windows.append(sequence[i:i + window_size])
    return windows


def create_af3_json(name, sequence):
    """Create AlphaFold3 JSON structure with a single model seed."""
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate AlphaFold3 JSON input files for MSA generation"
    )
    parser.add_argument("--fasta", type=str, required=True,
                        help="FASTA file containing protein sequence(s)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory for JSON files")
    parser.add_argument("--sliding", action="store_true",
                        help="Enable sliding window mode")
    parser.add_argument("-k", type=int, default=15,
                        help="Window size for sliding mode (default: 15)")
    parser.add_argument("--name", type=str,
                        help="Custom name for output (only for single sequence FASTA)")

    args = parser.parse_args()

    # Parse input FASTA
    sequences = parse_fasta(args.fasta)

    if len(sequences) == 0:
        print("Error: No sequences found in FASTA file")
        return

    # Check --name usage
    if args.name and len(sequences) > 1:
        print("Warning: --name option ignored for multi-sequence FASTA")
        args.name = None

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    total_files = 0

    for seq_id, seq in sequences.items():
        base_name = args.name if args.name else seq_id

        if args.sliding:
            windows = generate_sliding_windows(seq, args.k)

            for i, window_seq in enumerate(windows):
                json_name = f"{base_name}_{args.k}mer_window_{i:03d}"
                af3_json = create_af3_json(json_name, window_seq)

                filename = os.path.join(args.output_dir, f"{json_name}.json")
                with open(filename, "w") as f:
                    json.dump(af3_json, f, indent=4)

                total_files += 1

            print(f"Generated {len(windows)} sliding window files for {seq_id}")
        else:
            af3_json = create_af3_json(base_name, seq)

            filename = os.path.join(args.output_dir, f"{base_name}.json")
            with open(filename, "w") as f:
                json.dump(af3_json, f, indent=4)

            total_files += 1
            print(f"Generated file for {seq_id}")

    print(f"\nTotal files generated: {total_files}")
    print(f"Output directory: {args.output_dir}")
    if args.sliding:
        print(f"Window size: {args.k} AA")


if __name__ == "__main__":
    main()
