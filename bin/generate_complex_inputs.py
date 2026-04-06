#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files for specified peptides bound to antibodies.
Creates JSON files with peptide sequences (or sliding windows) and antibody heavy/light chains.
"""

import argparse
import json
import os
from pathlib import Path
from typing import Optional, Dict
from utils import (parse_fasta, create_af3_json, create_protein_sequence,
                   generate_sliding_windows, generate_model_seeds)


def load_msa_index(index_path: Optional[str]) -> Dict:
    """Load MSA index if provided."""
    if not index_path or not os.path.exists(index_path):
        return {}
    with open(index_path, 'r') as f:
        return json.load(f)


def get_msa_data(msa_index: Dict, seq_key: str) -> Optional[Dict]:
    """Get MSA data with case-insensitive lookup."""
    # Try exact match first
    if seq_key in msa_index:
        return msa_index[seq_key]

    # Try lowercase match
    seq_key_lower = seq_key.lower()
    for key, value in msa_index.items():
        if key.lower() == seq_key_lower:
            return value

    return None


def staged_seq_dir(msa_data: Dict, staged_msa_dir: Optional[str]) -> Optional[Path]:
    """Return the staged directory for this sequence, or None if no staged dir."""
    if not staged_msa_dir:
        return None
    return Path(staged_msa_dir) / Path(msa_data['msa_dir']).name


def msa_exists(msa_data: Dict, staged_msa_dir: Optional[str]) -> bool:
    """Check that the unpaired MSA exists, using the staged directory when available."""
    seq_dir = staged_seq_dir(msa_data, staged_msa_dir)
    if seq_dir is not None:
        return (seq_dir / 'unpaired_msa.a3m').exists()
    return os.path.exists(msa_data['unpaired_msa'])


def load_template_data(msa_data: Dict, staged_msa_dir: Optional[str] = None) -> list:
    """Load template information from metadata.json."""
    seq_dir = staged_seq_dir(msa_data, staged_msa_dir)
    if seq_dir is not None:
        metadata_path = seq_dir / 'metadata.json'
        templates_base = seq_dir / 'templates'
    else:
        metadata_path = Path(msa_data['msa_dir']) / 'metadata.json'
        templates_base = Path(msa_data['templates_dir'])

    if not metadata_path.exists():
        return []

    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read {metadata_path}: {e}")
        return []

    if not metadata.get('has_templates', False):
        return []

    templates = []
    for i, template_meta in enumerate(metadata.get('templates', [])):
        template_cif = templates_base / f'template_{i}.cif'
        if template_cif.exists():
            # mmcifPath must point to the published location for AF3 to resolve later
            published_cif = Path(msa_data['templates_dir']) / f'template_{i}.cif'
            templates.append({
                'mmcifPath': str(published_cif),
                'queryIndices': template_meta['queryIndices'],
                'templateIndices': template_meta['templateIndices']
            })

    return templates


def extract_antibody_chains(antibody_sequences):
    """Extract heavy and light chains from antibody sequences based on naming."""
    heavy_chain = None
    light_chain = None
    heavy_id = None
    light_id = None

    for seq_id, sequence in antibody_sequences.items():
        seq_id_lower = seq_id.lower()
        if 'heavy' in seq_id_lower or '_hc' in seq_id_lower or seq_id_lower.endswith('hc'):
            heavy_chain = sequence
            heavy_id = seq_id
        elif 'light' in seq_id_lower or '_lc' in seq_id_lower or seq_id_lower.endswith('lc'):
            light_chain = sequence
            light_id = seq_id

    if heavy_chain is None or light_chain is None:
        raise ValueError("Could not identify heavy and light chains in antibody sequences. "
                        "Ensure sequence names contain 'heavy'/'light' or 'HC'/'LC'.")

    return heavy_chain, light_chain, heavy_id, light_id


def main():
    parser = argparse.ArgumentParser(description="Generate AlphaFold3 JSON input files for peptide-antibody complexes")
    parser.add_argument("--peptide-fasta", type=str, required=True,
                        help="FASTA file containing peptide sequences")
    parser.add_argument("--antibody-fasta", type=str, required=True,
                        help="FASTA file containing antibody heavy and light chain sequences")
    parser.add_argument("--antibody-name", type=str, required=True,
                        help="Name for the antibody to use in output filenames")
    parser.add_argument("--output-dir", type=str, default="af3_inputs",
                        help="Output directory for JSON files (default: af3_inputs)")
    parser.add_argument("--sliding", action="store_true",
                        help="Enable sliding window mode")
    parser.add_argument("-k", type=int, default=15,
                        help="Window size for sliding mode (default: 15)")
    parser.add_argument("--num-seeds", type=int, default=100,
                        help="Number of model seeds to generate (default: 100)")
    parser.add_argument("--msa-index", type=str, default=None,
                        help="JSON file mapping sequence IDs to MSA paths (optional)")
    parser.add_argument("--msa-dir", type=str, default=None,
                        help="Staged MSA directory (Nextflow-guaranteed to exist before script runs)")

    args = parser.parse_args()

    # Generate consistent random seeds
    model_seeds = generate_model_seeds(num_seeds=args.num_seeds)

    # Load MSA index if provided
    msa_index = load_msa_index(args.msa_index)
    if msa_index:
        print(f"Loaded MSA index with {len(msa_index)} entries")

    staged_msa_dir = args.msa_dir if args.msa_dir else None

    # Parse input files
    peptide_sequences = parse_fasta(args.peptide_fasta)
    antibody_sequences = parse_fasta(args.antibody_fasta)

    # Extract heavy and light chains
    heavy_chain, light_chain, heavy_id, light_id = extract_antibody_chains(antibody_sequences)
    
    print(f"Found {len(peptide_sequences)} peptide sequences")
    print(f"Heavy chain length: {len(heavy_chain)} AA")
    print(f"Light chain length: {len(light_chain)} AA")
    
    # Create output directory
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
                window_name = f"{peptide_name}_{args.k}mer_window_{i:03d}"

                # Create peptide sequence with MSA paths
                peptide_entry = create_protein_sequence("A", window_seq)
                msa_data = get_msa_data(msa_index, window_name)
                if msa_data and msa_exists(msa_data, staged_msa_dir):
                    peptide_entry['protein']['unpairedMsaPath'] = msa_data['unpaired_msa']
                    peptide_entry['protein']['pairedMsaPath'] = msa_data['paired_msa']
                    templates = load_template_data(msa_data, staged_msa_dir)
                    peptide_entry['protein']['templates'] = templates if templates else []

                # Create heavy chain with MSA paths
                heavy_entry = create_protein_sequence("B", heavy_chain)
                heavy_key = heavy_id
                heavy_data = get_msa_data(msa_index, heavy_key)
                if heavy_data and msa_exists(heavy_data, staged_msa_dir):
                    heavy_entry['protein']['unpairedMsaPath'] = heavy_data['unpaired_msa']
                    heavy_entry['protein']['pairedMsaPath'] = heavy_data['paired_msa']
                    templates = load_template_data(heavy_data, staged_msa_dir)
                    heavy_entry['protein']['templates'] = templates if templates else []

                # Create light chain with MSA paths
                light_entry = create_protein_sequence("C", light_chain)
                light_key = light_id
                light_data = get_msa_data(msa_index, light_key)
                if light_data and msa_exists(light_data, staged_msa_dir):
                    light_entry['protein']['unpairedMsaPath'] = light_data['unpaired_msa']
                    light_entry['protein']['pairedMsaPath'] = light_data['paired_msa']
                    templates = load_template_data(light_data, staged_msa_dir)
                    light_entry['protein']['templates'] = templates if templates else []

                sequences = [peptide_entry, heavy_entry, light_entry]

                json_name = f"{peptide_name}_{args.antibody_name}_{args.k}mer_window_{i:02d}"
                af3_json = create_af3_json(json_name, sequences, model_seeds)
                
                filename = os.path.join(args.output_dir, f"{json_name}.json")
                with open(filename, "w") as f:
                    json.dump(af3_json, f, indent=4)
                
                total_files += 1
            
            print(f"Generated {len(windows)} sliding window files for {peptide_name}")
        
        else:
            # Default mode - full peptide sequence
            # Create peptide sequence with MSA paths
            peptide_entry = create_protein_sequence("A", peptide_seq)
            msa_data = get_msa_data(msa_index, peptide_name)
            if msa_data and msa_exists(msa_data, staged_msa_dir):
                peptide_entry['protein']['unpairedMsaPath'] = msa_data['unpaired_msa']
                peptide_entry['protein']['pairedMsaPath'] = msa_data['paired_msa']
                templates = load_template_data(msa_data, staged_msa_dir)
                peptide_entry['protein']['templates'] = templates if templates else []

            # Create heavy chain with MSA paths
            heavy_entry = create_protein_sequence("B", heavy_chain)
            heavy_key = heavy_id
            heavy_data = get_msa_data(msa_index, heavy_key)
            if heavy_data and msa_exists(heavy_data, staged_msa_dir):
                heavy_entry['protein']['unpairedMsaPath'] = heavy_data['unpaired_msa']
                heavy_entry['protein']['pairedMsaPath'] = heavy_data['paired_msa']
                templates = load_template_data(heavy_data, staged_msa_dir)
                heavy_entry['protein']['templates'] = templates if templates else []

            # Create light chain with MSA paths
            light_entry = create_protein_sequence("C", light_chain)
            light_key = light_id
            light_data = get_msa_data(msa_index, light_key)
            if light_data and msa_exists(light_data, staged_msa_dir):
                light_entry['protein']['unpairedMsaPath'] = light_data['unpaired_msa']
                light_entry['protein']['pairedMsaPath'] = light_data['paired_msa']
                templates = load_template_data(light_data, staged_msa_dir)
                light_entry['protein']['templates'] = templates if templates else []

            sequences = [peptide_entry, heavy_entry, light_entry]

            json_name = f"{peptide_name}_{args.antibody_name}"
            af3_json = create_af3_json(json_name, sequences, model_seeds)
            
            filename = os.path.join(args.output_dir, f"{json_name}.json")
            with open(filename, "w") as f:
                json.dump(af3_json, f, indent=4)
            
            total_files += 1
            print(f"Generated file for {peptide_name}")
    
    print(f"\nTotal files generated: {total_files}")
    print(f"Output directory: {args.output_dir}")
    if args.sliding:
        print(f"Window size: {args.k} AA")


if __name__ == "__main__":
    main()