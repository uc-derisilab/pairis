#!/usr/bin/env python3
"""Build index mapping sequence IDs to MSA file locations."""

import argparse
import json
from pathlib import Path


def build_msa_index(msa_root: str, output_json: str, published_base: str = None) -> None:
    """
    Scan MSA directory and build index.

    Expected structure:
    msa_root/
      ABLIM2_fragment_17_window_000/
        unpaired_msa.a3m
        paired_msa.a3m
        templates/
        metadata.json

    Args:
        msa_root: Root directory containing MSA subdirectories
        output_json: Output JSON file path
        published_base: Base path where MSAs will be published (if different from msa_root)
    """
    msa_path = Path(msa_root)
    index = {}

    for metadata_file in msa_path.rglob('metadata.json'):
        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)

            # Always use directory name - guaranteed unique (metadata 'id' is often just "A")
            seq_id = metadata_file.parent.name
            msa_dir_rel = metadata_file.parent.relative_to(msa_path)

            # Use published base path if provided, otherwise use absolute work dir path
            if published_base:
                msa_dir = Path(published_base) / msa_dir_rel
            else:
                msa_dir = metadata_file.parent.absolute()

            # Build index entry with paths
            index[seq_id] = {
                'msa_dir': str(msa_dir),
                'unpaired_msa': str(msa_dir / 'unpaired_msa.a3m'),
                'paired_msa': str(msa_dir / 'paired_msa.a3m'),
                'templates_dir': str(msa_dir / 'templates'),
                'has_templates': metadata.get('has_templates', False),
                'num_templates': metadata.get('num_templates', 0)
            }

        except Exception as e:
            print(f"Error processing {metadata_file}: {e}")

    with open(output_json, 'w') as f:
        json.dump(index, f, indent=2)

    print(f"Built MSA index with {len(index)} entries")


def main():
    parser = argparse.ArgumentParser(description="Build MSA index")
    parser.add_argument("--msa-dir", required=True, help="MSA root directory")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--published-base", help="Base path where MSAs will be published")
    args = parser.parse_args()

    build_msa_index(args.msa_dir, args.output, args.published_base)


if __name__ == "__main__":
    main()
