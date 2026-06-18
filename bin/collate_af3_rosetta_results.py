#!/usr/bin/env python3
"""
Collate AF3 and Rosetta results for PAIRIS antibody-peptide binding analysis.

This script combines iPTM scores from AlphaFold3 predictions with binding energies
from Rosetta docking simulations for antibody-peptide complexes.
"""

import argparse
import json
import pandas as pd
from pathlib import Path
import re


def extract_peptide_sequence(af3_input_json_path):
    """Extract peptide (chain A) sequence from a folding-backend input JSON.

    Handles both input shapes:
    - AlphaFold3 (nested): {"sequences": [{"protein": {"id": "A", "sequence": ...}}]}
    - ESMFold2 (flat):      {"sequences": [{"type": "protein", "id": "A", "sequence": ...}]}

    Args:
        af3_input_json_path: Path to the folding-backend input JSON file

    Returns:
        Peptide sequence string, or None if not found
    """
    try:
        with open(af3_input_json_path, 'r') as f:
            data = json.load(f)

        # Find the protein sequence with id "A" (the peptide) in either shape.
        for seq_entry in data.get('sequences', []):
            # AF3 nested shape
            if 'protein' in seq_entry:
                if seq_entry['protein'].get('id') == 'A':
                    return seq_entry['protein'].get('sequence')
            # ESMFold2 flat shape
            elif seq_entry.get('type') == 'protein' and seq_entry.get('id') == 'A':
                return seq_entry.get('sequence')
        return None
    except Exception as e:
        print(f"Warning: Error reading {af3_input_json_path}: {e}")
        return None


def extract_iptm(af3_output_json_path):
    """Extract iPTM score from AF3 output summary confidences JSON.

    Args:
        af3_output_json_path: Path to summary_confidences.json file

    Returns:
        iPTM score (float), or None if not found
    """
    try:
        with open(af3_output_json_path, 'r') as f:
            data = json.load(f)

        # Get first element of chain_iptm
        chain_iptm = data.get('chain_iptm', [])
        if chain_iptm:
            return chain_iptm[0]
        return None
    except Exception as e:
        print(f"Warning: Error reading {af3_output_json_path}: {e}")
        return None


def parse_pairis_filename(filename):
    """Parse peptide name, BCR, and window from PAIRIS filename.

    Args:
        filename: String like "BCAS3_epitope_mAb1_HA1C_15mer_window_00" or
                  "BCAS3_epitope_rosetta.csv" or
                  "BCAS3_epitope_mAb1_HA1C_15mer_window_00_bcas3_epitope_mab1_ha1c_15mer_window_00_model"

    Returns:
        Tuple of (peptide_name, bcr_name, window) where window may be None
    """
    # Pattern: {PEPTIDE}_{BCR}_{KMER}mer_window_{N}[_anything_else]
    # This handles both simple names and complex Rosetta output names
    pattern = r'^(.+?)_([^_]+_[^_]+)_\d+mer_window_(\d+)'
    match = re.search(pattern, filename)

    if match:
        peptide_name = match.group(1)
        bcr_name = match.group(2)
        window = match.group(3)
        return peptide_name, bcr_name, window

    # Try pattern without window (for Rosetta CSV filenames)
    pattern_no_window = r'^(.+)_rosetta$'
    match = re.match(pattern_no_window, filename)
    if match:
        peptide_name = match.group(1)
        return peptide_name, None, None

    return None, None, None


def find_af3_output_json(af3_output_dir, kmer_size, bcr, peptide, window):
    """Find the summary_confidences.json file for a specific run.

    Args:
        af3_output_dir: Base AF3 output directory
        kmer_size: Kmer size (e.g., "15mers" - plural from directory structure)
        bcr: BCR identifier (e.g., "mAb1_HA1C")
        peptide: Peptide name (e.g., "BCAS3_epitope")
        window: Window number (as string, e.g., "00")

    Returns:
        Path to summary_confidences.json or None if not found
    """
    # Convert kmer_size from directory format "15mers" to filename format "15mer"
    kmer_num = kmer_size.replace('mers', 'mer')

    # Construct the directory name - PAIRIS uses uppercase for outer, lowercase for inner
    outer_dir_name = f"{peptide}_{bcr}_{kmer_num}_window_{window}"
    inner_dir_name = outer_dir_name.lower()

    # The path structure in PAIRIS is: complexes/{PEPTIDE}_{BCR}_{KMER}mer_window_{N}/{inner}/{files}
    base_path = Path(af3_output_dir) / outer_dir_name / inner_dir_name

    if not base_path.exists():
        return None

    json_path = base_path / f"{inner_dir_name}_summary_confidences.json"
    if json_path.exists():
        return json_path

    return None


def collate_bcr_data(af3_input_dir, af3_output_dir, rosetta_dir, bcr, kmer_size):
    """Collate all data for a single BCR and kmer size.

    Args:
        af3_input_dir: Path to AF3 inputs directory
        af3_output_dir: Path to AF3 outputs directory
        rosetta_dir: Path to Rosetta results directory
        bcr: BCR identifier
        kmer_size: Kmer size (e.g., "15mers")

    Returns:
        List of dictionaries containing collated data
    """
    results = []

    # Get all Rosetta CSV files for this BCR
    rosetta_bcr_dir = Path(rosetta_dir) / kmer_size / bcr

    if not rosetta_bcr_dir.exists():
        print(f"Warning: Rosetta directory not found: {rosetta_bcr_dir}")
        return results

    rosetta_csv_files = list(rosetta_bcr_dir.glob("*_rosetta.csv"))

    # Full-length run (window_sizes: null): no window, one structure per BCR-peptide
    if kmer_size == 'full':
        for csv_file in rosetta_csv_files:
            peptide_name, _, _ = parse_pairis_filename(csv_file.stem)
            if not peptide_name:
                print(f"Warning: Could not parse filename: {csv_file.name}")
                continue
            try:
                rosetta_df = pd.read_csv(csv_file)
            except Exception as e:
                print(f"Error reading {csv_file}: {e}")
                continue
            for _, row in rosetta_df.iterrows():
                binding_energy = row.get('binding_energy', None)
                outer_dir = f"{peptide_name}_{bcr}"
                inner_dir = outer_dir.lower()
                af3_input_json = Path(af3_input_dir) / f"{outer_dir}.json"
                peptide_seq = extract_peptide_sequence(af3_input_json) if af3_input_json.exists() else None
                base_path = Path(af3_output_dir) / outer_dir / inner_dir
                af3_output_json = base_path / f"{inner_dir}_summary_confidences.json"
                iptm = extract_iptm(af3_output_json) if af3_output_json.exists() else None
                results.append({
                    'peptide_name': peptide_name,
                    'window': None,
                    'peptide_sequence': peptide_seq,
                    'kmer_size': None,
                    'iptm': iptm,
                    'binding_energy': binding_energy
                })
        return results

    for csv_file in rosetta_csv_files:
        # Parse peptide name from filename
        peptide_name, _, _ = parse_pairis_filename(csv_file.stem)

        if not peptide_name:
            print(f"Warning: Could not parse filename: {csv_file.name}")
            continue

        # Load Rosetta CSV
        try:
            rosetta_df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

        # Process each window in the Rosetta CSV
        for _, row in rosetta_df.iterrows():
            # Parse window number from filename column
            _, _, window = parse_pairis_filename(row['filename'])

            if window is None:
                print(f"Warning: Could not parse window from: {row['filename']}")
                continue

            # Get binding energy
            binding_energy = row.get('binding_energy', None)

            # Convert kmer_size from directory format "15mers" to filename format "15mer"
            kmer_num = kmer_size.replace('mers', 'mer')

            # Find AF3 input JSON to get peptide sequence
            af3_input_json = Path(af3_input_dir) / f"{peptide_name}_{bcr}_{kmer_num}_window_{window}.json"
            peptide_seq = extract_peptide_sequence(af3_input_json) if af3_input_json.exists() else None

            # Find AF3 output JSON to get iPTM
            af3_output_json = find_af3_output_json(af3_output_dir, kmer_size, bcr, peptide_name, window)
            iptm = extract_iptm(af3_output_json) if af3_output_json else None

            if af3_output_json is None:
                print(f"Warning: Could not find summary_confidences.json for {peptide_name}_{bcr}_{kmer_num}_window_{window}")

            # Add to results
            results.append({
                'peptide_name': peptide_name,
                'window': int(window),
                'peptide_sequence': peptide_seq,
                'kmer_size': int(kmer_size.replace('mers', '')),
                'iptm': iptm,
                'binding_energy': binding_energy
            })

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Collate AF3 and Rosetta results for PAIRIS antibody-peptide binding analysis'
    )
    parser.add_argument('--af3-input-dir', required=True,
                        help='Directory containing AF3 input JSONs')
    parser.add_argument('--af3-output-dir', required=True,
                        help='Directory containing AF3 output results')
    parser.add_argument('--rosetta-dir', required=True,
                        help='Directory containing Rosetta results')
    parser.add_argument('--output-dir', required=True,
                        help='Directory to write output CSV files')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get list of BCRs from Rosetta directory
    rosetta_base = Path(args.rosetta_dir)
    bcrs = set()

    # Find all kmer size directories
    kmer_sizes = []
    for item in rosetta_base.iterdir():
        if item.is_dir() and ('mers' in item.name or item.name == 'full'):
            kmer_sizes.append(item.name)
            # Get BCRs from this kmer directory
            for bcr_dir in item.iterdir():
                if bcr_dir.is_dir():
                    bcrs.add(bcr_dir.name)

    bcrs = sorted(bcrs)
    kmer_sizes = sorted(kmer_sizes)

    print(f"Found kmer sizes: {', '.join(kmer_sizes)}")
    print(f"Found BCRs: {', '.join(bcrs)}")

    # Process each BCR
    for bcr in bcrs:
        print(f"\nProcessing {bcr}...")
        all_results = []

        # Process all kmer sizes
        for kmer_size in kmer_sizes:
            print(f"  Processing {kmer_size}...")
            results = collate_bcr_data(
                af3_input_dir=args.af3_input_dir,
                af3_output_dir=args.af3_output_dir,
                rosetta_dir=args.rosetta_dir,
                bcr=bcr,
                kmer_size=kmer_size
            )
            all_results.extend(results)
            print(f"    Found {len(results)} entries")

        # Convert to DataFrame and save
        if all_results:
            df = pd.DataFrame(all_results)
            # Reorder columns
            df = df[['peptide_name', 'window', 'peptide_sequence', 'kmer_size', 'iptm', 'binding_energy']]
            # Sort by peptide name, window
            df = df.sort_values(['peptide_name', 'window'])

            output_file = output_dir / f"{bcr}.csv"
            df.to_csv(output_file, index=False)
            print(f"  Wrote {len(df)} rows to {output_file}")
        else:
            print(f"  No data found for {bcr}")

    print("\nDone!")


if __name__ == '__main__':
    main()
