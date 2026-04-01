#!/usr/bin/env python3
"""
Shared utility functions for AlphaFold3 input generation scripts.
"""

import random
from typing import List, Dict, Any


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """Parse FASTA file and return sequences as dictionary."""
    sequences = {}
    current_id = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]  # Remove '>'
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line
    
    return sequences


def create_af3_json(name: str, sequences: List[Dict[str, Any]], seeds: List[int]) -> Dict[str, Any]:
    """Create AlphaFold3 JSON structure."""
    af3_json = {
        "name": name,
        "modelSeeds": seeds,
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 2  # Version 2 supports MSA paths
    }

    return af3_json


def create_protein_sequence(seq_id: str, sequence: str) -> Dict[str, Any]:
    """Create protein sequence entry for AF3 JSON."""
    protein_entry = {
        "protein": {
            "id": seq_id,
            "sequence": sequence
        }
    }
    
    return protein_entry


def generate_sliding_windows(sequence: str, window_size: int = 15) -> List[str]:
    """Generate sliding windows of specified size across sequence."""
    windows = []
    for i in range(len(sequence) - window_size + 1):
        windows.append(sequence[i:i + window_size])
    return windows


def generate_model_seeds(seed: int = 42, num_seeds: int = 100) -> List[int]:
    """Generate consistent random seeds for AlphaFold3 modeling."""
    random.seed(seed)
    return [random.randint(1, 10000) for _ in range(num_seeds)]