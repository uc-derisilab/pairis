#!/usr/bin/env python3
"""
Generate ESMFold2 JobSpec JSON input files for peptides bound to antibodies.

Mirrors generate_complex_inputs.py (the AlphaFold3 generator) so that output
filenames stay byte-identical across folding backends. Produces one JSON per
(peptide x antibody x window) combination with peptide/heavy/light chains.

The ESMFold2 schema has no template support and only the "full" model accepts
per-chain MSAs; the "fast" model rejects any `msas` key.
"""

import argparse
import json
import os
from utils import parse_fasta, generate_sliding_windows
# Reuse the AF3 generator helpers so naming/MSA logic stays identical.
from generate_complex_inputs import (
    extract_antibody_chains,
    load_msa_index,
    get_msa_data,
    msa_exists,
)


def create_protein_sequence(seq_id: str, sequence: str) -> dict:
    """Create an ESMFold2 protein sequence entry."""
    return {
        "type": "protein",
        "id": seq_id,
        "sequence": sequence,
    }


def make_msa_entry(chain_id: str, msa_data: dict) -> dict:
    """Create an ESMFold2 MSA entry pointing at the unpaired MSA a3m."""
    return {
        "chain_id": chain_id,
        "a3m_path": msa_data["unpaired_msa"],
    }


def create_jobspec(name: str, sequences: list, config: dict,
                   msas: list = None) -> dict:
    """Build an ESMFold2 JobSpec dict. Omit `msas` entirely when empty."""
    spec = {
        "schema_version": 1,
        "name": name,
        "sequences": sequences,
        "config": config,
    }
    if msas:
        spec["msas"] = msas
    return spec


def build_msas(use_msas: bool, chain_lookup, msa_index: dict,
               staged_msa_dir) -> list:
    """Build the msas list for a complex.

    chain_lookup: list of (chain_id, lookup_key) in chain order.
    Returns [] when use_msas is False or no chain has an existing MSA.
    """
    if not use_msas:
        return []
    msas = []
    for chain_id, lookup_key in chain_lookup:
        msa_data = get_msa_data(msa_index, lookup_key)
        if msa_data and msa_exists(msa_data, staged_msa_dir):
            msas.append(make_msa_entry(chain_id, msa_data))
    return msas


def main():
    parser = argparse.ArgumentParser(
        description="Generate ESMFold2 JobSpec JSON input files for "
                    "peptide-antibody complexes")
    parser.add_argument("--peptide-fasta", type=str, required=True,
                        help="FASTA file containing peptide sequences")
    parser.add_argument("--antibody-fasta", type=str, required=True,
                        help="FASTA file containing antibody heavy and light "
                             "chain sequences")
    parser.add_argument("--antibody-name", type=str, required=True,
                        help="Name for the antibody to use in output filenames")
    parser.add_argument("--output-dir", type=str, default=".",
                        help="Output directory for JSON files (default: .)")
    parser.add_argument("--sliding", action="store_true",
                        help="Enable sliding window mode")
    parser.add_argument("-k", type=int, default=15,
                        help="Window size for sliding mode (default: 15)")
    parser.add_argument("--model", choices=["fast", "full"], default="fast",
                        help="ESMFold2 model variant (default: fast). "
                             "Only 'full' accepts MSAs.")
    parser.add_argument("--msa-index", type=str, default=None,
                        help="JSON file mapping sequence IDs to MSA paths "
                             "(optional; only used with --model full)")
    parser.add_argument("--msa-dir", type=str, default=None,
                        help="Staged MSA directory (Nextflow-guaranteed to "
                             "exist before script runs)")
    parser.add_argument("--num-loops", type=int, default=3,
                        help="ESMFold2 num_loops (default: 3)")
    parser.add_argument("--num-sampling-steps", type=int, default=50,
                        help="ESMFold2 num_sampling_steps (default: 50)")

    args = parser.parse_args()

    # The fast model rejects MSAs; only resolve them for the full model.
    use_msas = args.model == "full"

    msa_index = {}
    if use_msas:
        msa_index = load_msa_index(args.msa_index)
        if msa_index:
            print(f"Loaded MSA index with {len(msa_index)} entries")
    staged_msa_dir = args.msa_dir if args.msa_dir else None

    # Parse input files
    peptide_sequences = parse_fasta(args.peptide_fasta)
    antibody_sequences = parse_fasta(args.antibody_fasta)

    # Extract heavy and light chains
    heavy_chain, light_chain, heavy_id, light_id = extract_antibody_chains(
        antibody_sequences)

    print(f"Found {len(peptide_sequences)} peptide sequences")
    print(f"Heavy chain length: {len(heavy_chain)} AA")
    print(f"Light chain length: {len(light_chain)} AA")

    config = {
        "num_loops": args.num_loops,
        "num_sampling_steps": args.num_sampling_steps,
        "seed": 0,
    }

    os.makedirs(args.output_dir, exist_ok=True)

    total_files = 0

    for peptide_name, peptide_seq in peptide_sequences.items():
        if args.sliding:
            # Sliding window mode
            if len(peptide_seq) < args.k:
                # Use full sequence if shorter than window size
                windows = [peptide_seq]
            else:
                windows = generate_sliding_windows(peptide_seq, args.k)

            for i, window_seq in enumerate(windows):
                # 3-digit form is the MSA-lookup key for the peptide chain.
                window_name = f"{peptide_name}_{args.k}mer_window_{i:03d}"

                sequences = [
                    create_protein_sequence("A", window_seq),
                    create_protein_sequence("B", heavy_chain),
                    create_protein_sequence("C", light_chain),
                ]

                msas = build_msas(use_msas, [
                    ("A", window_name),
                    ("B", heavy_id),
                    ("C", light_id),
                ], msa_index, staged_msa_dir)

                # 2-digit form is the output name AND file stem.
                json_name = (f"{peptide_name}_{args.antibody_name}_"
                             f"{args.k}mer_window_{i:02d}")
                spec = create_jobspec(json_name, sequences, config, msas)

                filename = os.path.join(args.output_dir, f"{json_name}.json")
                with open(filename, "w") as f:
                    json.dump(spec, f, indent=4)

                total_files += 1

            print(f"Generated {len(windows)} sliding window files "
                  f"for {peptide_name}")

        else:
            # Default mode - full peptide sequence
            sequences = [
                create_protein_sequence("A", peptide_seq),
                create_protein_sequence("B", heavy_chain),
                create_protein_sequence("C", light_chain),
            ]

            msas = build_msas(use_msas, [
                ("A", peptide_name),
                ("B", heavy_id),
                ("C", light_id),
            ], msa_index, staged_msa_dir)

            json_name = f"{peptide_name}_{args.antibody_name}"
            spec = create_jobspec(json_name, sequences, config, msas)

            filename = os.path.join(args.output_dir, f"{json_name}.json")
            with open(filename, "w") as f:
                json.dump(spec, f, indent=4)

            total_files += 1
            print(f"Generated file for {peptide_name}")

    print(f"\nTotal files generated: {total_files}")
    print(f"Output directory: {args.output_dir}")
    if args.sliding:
        print(f"Window size: {args.k} AA")


if __name__ == "__main__":
    main()
