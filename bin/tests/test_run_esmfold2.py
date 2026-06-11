#!/usr/bin/env python3
"""Tests for run_esmfold2.py pure helpers (no GPU / torch / esm / transformers).

These tests only exercise the GPU-free helper functions. The model-loading,
input-building, and folding functions import torch/esm/transformers lazily
inside their bodies, so importing this module must NOT require those packages.
"""

import importlib
import sys
from pathlib import Path

import pytest

# Make bin/ importable (this file lives in bin/tests/)
BIN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BIN_DIR))

import run_esmfold2  # noqa: E402


def test_import_has_no_heavy_imports():
    """Importing run_esmfold2 must not pull in torch/esm/transformers.

    The module imports cleanly in an env without GPU libs, and the pure
    helper is callable.
    """
    mod = importlib.import_module("run_esmfold2")
    assert callable(mod.compute_chain_iptm)
    # None of the heavy GPU libraries should have been imported as a side
    # effect of importing the module.
    assert "torch" not in sys.modules
    assert "esm" not in sys.modules
    assert "transformers" not in sys.modules


def test_compute_chain_iptm_3x3():
    # chain_iptm[i] = mean over j != i of m[i][j]
    m = [
        [0.0, 0.2, 0.4],   # row 0: mean(0.2, 0.4) = 0.3
        [0.6, 0.0, 0.8],   # row 1: mean(0.6, 0.8) = 0.7
        [0.1, 0.5, 0.0],   # row 2: mean(0.1, 0.5) = 0.3
    ]
    result = run_esmfold2.compute_chain_iptm(m)
    assert result == pytest.approx([0.3, 0.7, 0.3])


def test_compute_chain_iptm_single_chain():
    # n_chains == 1 -> no off-diagonal -> [0.0]
    assert run_esmfold2.compute_chain_iptm([[0.95]]) == [0.0]


def test_compute_chain_iptm_nested_lists_accepted():
    # Plain nested lists must work (don't assume numpy in the signature).
    m = [[0.0, 1.0], [1.0, 0.0]]
    assert run_esmfold2.compute_chain_iptm(m) == pytest.approx([1.0, 1.0])


def test_output_paths_lowercases_inner():
    outer = "BCAS3_epitope_mAb1_HA1C_15mer_window_00"
    inner = "bcas3_epitope_mab1_ha1c_15mer_window_00"
    cif, summary, all_seeds = run_esmfold2.output_paths("/some/out", outer)

    assert cif == Path("/some/out") / inner / f"{inner}_model.cif"
    assert summary == Path("/some/out") / inner / f"{inner}_summary_confidences.json"
    assert all_seeds == Path("/some/out") / inner / f"{inner}_all_seeds.json"


def test_output_paths_nesting_and_filenames():
    outer = "Pep1_MAB1"
    inner = "pep1_mab1"
    cif, summary, all_seeds = run_esmfold2.output_paths(".", outer)

    # Nested under <output_dir>/<inner>/
    assert cif.parent.name == inner
    assert summary.parent.name == inner
    assert all_seeds.parent.name == inner
    # Exact filenames
    assert cif.name == f"{inner}_model.cif"
    assert summary.name == f"{inner}_summary_confidences.json"
    assert all_seeds.name == f"{inner}_all_seeds.json"


def test_build_summary_confidences_preserves_order():
    chain_iptm = [0.3, 0.7, 0.5]  # peptide first
    d = run_esmfold2.build_summary_confidences(
        chain_iptm=chain_iptm, ptm=0.82, iptm=0.61, best_seed=4
    )
    # chain_iptm preserved in order (peptide first)
    assert d["chain_iptm"] == [0.3, 0.7, 0.5]
    assert d["ptm"] == 0.82
    assert d["iptm"] == 0.61
    # ranking_score is the peptide interface value (chain_iptm[0])
    assert d["ranking_score"] == 0.3
    assert d["best_seed"] == 4


def test_pick_best_by_chain_iptm0():
    records = [
        {"seed": 0, "chain_iptm": [0.3, 0.7, 0.5], "ptm": 0.9, "iptm": 0.5},
        {"seed": 1, "chain_iptm": [0.6, 0.2, 0.4], "ptm": 0.4, "iptm": 0.5},
        {"seed": 2, "chain_iptm": [0.5, 0.5, 0.5], "ptm": 0.8, "iptm": 0.5},
    ]
    best = run_esmfold2.pick_best(records)
    assert best["seed"] == 1  # highest chain_iptm[0] = 0.6


def test_pick_best_tie_break_by_ptm():
    records = [
        {"seed": 0, "chain_iptm": [0.5, 0.7, 0.5], "ptm": 0.70, "iptm": 0.5},
        {"seed": 1, "chain_iptm": [0.5, 0.2, 0.4], "ptm": 0.90, "iptm": 0.5},
        {"seed": 2, "chain_iptm": [0.5, 0.5, 0.5], "ptm": 0.80, "iptm": 0.5},
    ]
    # All tie on chain_iptm[0] = 0.5 -> highest ptm wins (seed 1)
    best = run_esmfold2.pick_best(records)
    assert best["seed"] == 1


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
