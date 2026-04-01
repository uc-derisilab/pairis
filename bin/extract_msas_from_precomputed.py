#!/usr/bin/env python3
"""
Extract MSAs and templates from pre-computed AlphaFold3 JSON files.
Saves them as separate A3M and mmCIF files for reuse.
"""

import argparse
import json
import os
from pathlib import Path
from multiprocessing import Pool
from typing import Dict, Tuple


def extract_msa_data(input_json_path: str, output_dir: str) -> Tuple[str, bool]:
    """
    Extract MSAs and templates from a single JSON file.
    Returns (file_path, success).
    """
    try:
        with open(input_json_path, 'r') as f:
            data = json.load(f)

        # Get the first (and typically only) sequence
        if not data.get('sequences') or not data['sequences']:
            return (input_json_path, False)

        sequence_data = data['sequences'][0]
        protein_data = sequence_data.get('protein', {})

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Save unpaired MSA
        unpaired_msa = protein_data.get('unpairedMsa', '')
        if unpaired_msa:
            with open(os.path.join(output_dir, 'unpaired_msa.a3m'), 'w') as f:
                f.write(unpaired_msa)

        # Save paired MSA
        paired_msa = protein_data.get('pairedMsa', '')
        if paired_msa:
            with open(os.path.join(output_dir, 'paired_msa.a3m'), 'w') as f:
                f.write(paired_msa)

        # Save templates if they exist
        templates = protein_data.get('templates', [])
        template_metadata = []
        if templates:
            template_dir = os.path.join(output_dir, 'templates')
            os.makedirs(template_dir, exist_ok=True)

            for i, template in enumerate(templates):
                mmcif_data = template.get('mmcif', '')
                if mmcif_data:
                    with open(os.path.join(template_dir, f'template_{i}.cif'), 'w') as f:
                        f.write(mmcif_data)

                # Save template metadata (queryIndices and templateIndices)
                template_metadata.append({
                    'queryIndices': template.get('queryIndices', []),
                    'templateIndices': template.get('templateIndices', [])
                })

        # Save metadata
        metadata = {
            'sequence': protein_data.get('sequence', ''),
            'id': protein_data.get('id', ''),
            'name': data.get('name', ''),
            'has_templates': len(templates) > 0,
            'num_templates': len(templates),
            'templates': template_metadata
        }
        with open(os.path.join(output_dir, 'metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)

        return (input_json_path, True)

    except Exception as e:
        print(f"Error processing {input_json_path}: {e}")
        return (input_json_path, False)


def find_json_files(input_dir: str, output_dir: str):
    """
    Recursively find all *_data.json files and create corresponding output paths.
    """
    tasks = []
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    for json_file in input_path.rglob('*_data.json'):
        # Use the JSON filename stem to create unique output directory
        # This works regardless of whether files are symlinked flat or nested
        # Example: bcas3_epitope_window_000_data.json -> bcas3_epitope_window_000/
        seq_id = json_file.stem.replace('_data', '')
        output_subdir = output_path / seq_id

        tasks.append((str(json_file), str(output_subdir)))

    return tasks


def validate_extraction(output_dir: str) -> None:
    """Validate that each MSA directory is unique."""
    output_path = Path(output_dir)
    msa_dirs = set()
    duplicates = []

    for msa_file in output_path.rglob('unpaired_msa.a3m'):
        dir_name = msa_file.parent.name
        if dir_name in msa_dirs:
            duplicates.append(dir_name)
        msa_dirs.add(dir_name)

    if duplicates:
        print(f"WARNING: {len(duplicates)} duplicate MSA directories detected")
        for dup in duplicates[:5]:  # Show first 5
            print(f"  - {dup}")
    else:
        print(f"Validation passed: {len(msa_dirs)} unique MSA directories created")


def worker(args):
    """Worker function for multiprocessing."""
    return extract_msa_data(*args)


def main():
    parser = argparse.ArgumentParser(description="Extract MSAs and templates from pre-computed AF3 JSON files")
    parser.add_argument("--input-dir", type=str, required=True,
                        help="Directory containing AF3 output JSON files")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory for extracted MSA files")
    parser.add_argument("--num-workers", type=int, default=8,
                        help="Number of parallel workers (default: 8)")

    args = parser.parse_args()

    # Collect all tasks
    print("Collecting JSON files to process...")
    all_tasks = find_json_files(args.input_dir, args.output_dir)

    print(f"Found {len(all_tasks)} JSON files to process")

    # Process in parallel
    print(f"\nProcessing with {args.num_workers} workers...")
    with Pool(args.num_workers) as pool:
        results = pool.map(worker, all_tasks)

    # Summary
    successful = sum(1 for _, success in results if success)
    failed = len(results) - successful

    print(f"\nComplete!")
    print(f"Successfully processed: {successful}/{len(results)}")
    if failed > 0:
        print(f"Failed: {failed}")
        print("\nFailed files:")
        for path, success in results:
            if not success:
                print(f"  {path}")

    # Validate extraction results
    print(f"\nValidating extraction...")
    validate_extraction(args.output_dir)

    print(f"\nOutput saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
