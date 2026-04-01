#!/usr/bin/env python3
"""
Group structure directories by (kmer_size, BCR, peptide).
"""
import os
import re
import sys
from pathlib import Path

def main():
    if len(sys.argv) != 4:
        print("Usage: group_structures_by_bcr.py <complexes_dir> <bcr_dir> <window_size>")
        sys.exit(1)

    complexes_dir = sys.argv[1]
    bcr_dir = sys.argv[2]
    # window_size arg kept for API compatibility but no longer used in logic

    # Get BCR names from input directory
    bcr_names = []
    if os.path.exists(bcr_dir):
        for f in os.listdir(bcr_dir):
            if f.endswith('.fasta'):
                bcr_names.append(f.replace('.fasta', ''))

    bcr_set = set(bcr_names)  # O(1) lookups

    # Find all structure directories
    dirs = [d for d in os.listdir(complexes_dir)
            if os.path.isdir(os.path.join(complexes_dir, d))]

    # One compiled regex to strip the kmer+window suffix
    suffix_re = re.compile(r'_(\d+)mer_window_\d+$')

    groups = set()

    for dir_name in dirs:
        m = suffix_re.search(dir_name)
        if not m:
            continue
        kmer = m.group(1)
        prefix = dir_name[:m.start()]   # e.g. "antiHu_pep_108_HC_100_LC"

        # Try each split point left-to-right until a BCR match is found
        parts = prefix.split('_')
        for i in range(1, len(parts)):
            candidate_bcr = '_'.join(parts[i:])
            if candidate_bcr in bcr_set:
                peptide = '_'.join(parts[:i])
                groups.add((f"{kmer}mers", candidate_bcr, peptide))
                break

    # Write groups to file (one per line: kmer_size,bcr,peptide,complexes_dir)
    with open('groups.txt', 'w') as f:
        for kmer_size, bcr, peptide in sorted(groups):
            f.write(f"{kmer_size},{bcr},{peptide},{complexes_dir}\n")

if __name__ == '__main__':
    main()
