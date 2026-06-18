#!/usr/bin/env python3
"""Tests for collate_af3_rosetta_results.extract_peptide_sequence.

Covers both folding-backend input shapes: AlphaFold3 (nested under 'protein')
and ESMFold2 (flat protein entry with type/id at the top level).
"""

import json
import sys
from pathlib import Path

# Make bin/ importable (this file lives in bin/tests/)
BIN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BIN_DIR))

from collate_af3_rosetta_results import extract_peptide_sequence  # noqa: E402

PEPTIDE_SEQ = "EQKLISEEDLN"
HEAVY_SEQ = "QVQLVESGGGLVQPGGSLRLSCAAS"
LIGHT_SEQ = "DIQMTQSPSSLSASVGDRVTITCRAS"


def _write(tmp_path, data):
    p = tmp_path / "input.json"
    p.write_text(json.dumps(data))
    return p


def test_af3_nested_shape(tmp_path):
    data = {
        "sequences": [
            {"protein": {"id": "A", "sequence": PEPTIDE_SEQ}},
            {"protein": {"id": "B", "sequence": HEAVY_SEQ}},
            {"protein": {"id": "C", "sequence": LIGHT_SEQ}},
        ]
    }
    assert extract_peptide_sequence(_write(tmp_path, data)) == PEPTIDE_SEQ


def test_esmfold2_flat_shape(tmp_path):
    data = {
        "schema_version": 1,
        "name": "pep1_mab1",
        "sequences": [
            {"type": "protein", "id": "A", "sequence": PEPTIDE_SEQ},
            {"type": "protein", "id": "B", "sequence": HEAVY_SEQ},
            {"type": "protein", "id": "C", "sequence": LIGHT_SEQ},
        ],
    }
    assert extract_peptide_sequence(_write(tmp_path, data)) == PEPTIDE_SEQ


def test_no_chain_a_returns_none(tmp_path):
    data = {
        "sequences": [
            {"type": "protein", "id": "B", "sequence": HEAVY_SEQ},
        ]
    }
    assert extract_peptide_sequence(_write(tmp_path, data)) is None


def test_missing_file_returns_none(tmp_path):
    assert extract_peptide_sequence(tmp_path / "does_not_exist.json") is None
