#!/usr/bin/env python3
"""Tests for generate_esmfold2_inputs.py (ESMFold2 JobSpec generation)."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Make bin/ importable (this file lives in bin/tests/)
BIN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BIN_DIR))

from utils import generate_sliding_windows  # noqa: E402

SCRIPT = BIN_DIR / "generate_esmfold2_inputs.py"

HEAVY_SEQ = "QVQLVESGGGLVQPGGSLRLSCAAS"
LIGHT_SEQ = "DIQMTQSPSSLSASVGDRVTITCRAS"
PEPTIDE_SEQ = "ACDEFGHIKLMNPQRSTVWY"  # 20 aa


def write_fasta(path, records):
    """records: list of (id, seq)."""
    with open(path, "w") as f:
        for rid, seq in records:
            f.write(f">{rid}\n{seq}\n")
    return str(path)


def run_script(args):
    """Run the generator as a subprocess, return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
    )


def make_peptide_fasta(tmp_path, name="pep1", seq=PEPTIDE_SEQ):
    return write_fasta(tmp_path / "peptide.fasta", [(name, seq)])


def make_antibody_fasta(tmp_path, heavy_id="mab1_heavy", light_id="mab1_light"):
    return write_fasta(
        tmp_path / "antibody.fasta",
        [(heavy_id, HEAVY_SEQ), (light_id, LIGHT_SEQ)],
    )


def test_full_length_fast_no_msas(tmp_path):
    pep = make_peptide_fasta(tmp_path, name="pep1")
    ab = make_antibody_fasta(tmp_path)
    out = tmp_path / "out"

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--model", "fast",
    ])
    assert res.returncode == 0, res.stderr

    expected = out / "pep1_mab1.json"
    assert expected.exists(), f"missing {expected}; got {list(out.iterdir())}"

    spec = json.loads(expected.read_text())
    assert spec["schema_version"] == 1
    assert spec["name"] == "pep1_mab1"
    assert "msas" not in spec, "fast model must not emit msas"

    seqs = spec["sequences"]
    assert len(seqs) == 3
    assert [s["id"] for s in seqs] == ["A", "B", "C"]
    assert [s["type"] for s in seqs] == ["protein", "protein", "protein"]
    assert seqs[0]["sequence"] == PEPTIDE_SEQ
    assert seqs[1]["sequence"] == HEAVY_SEQ
    assert seqs[2]["sequence"] == LIGHT_SEQ

    cfg = spec["config"]
    assert cfg["num_loops"] == 3
    assert cfg["num_sampling_steps"] == 50
    assert cfg["seed"] == 0


def test_config_values_passed_through(tmp_path):
    pep = make_peptide_fasta(tmp_path, name="pep1")
    ab = make_antibody_fasta(tmp_path)
    out = tmp_path / "out"

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--model", "fast",
        "--num-loops", "7",
        "--num-sampling-steps", "13",
    ])
    assert res.returncode == 0, res.stderr
    spec = json.loads((out / "pep1_mab1.json").read_text())
    assert spec["config"]["num_loops"] == 7
    assert spec["config"]["num_sampling_steps"] == 13
    assert spec["config"]["seed"] == 0


def test_sliding_longer_than_k(tmp_path):
    k = 15
    pep = make_peptide_fasta(tmp_path, name="pep1", seq=PEPTIDE_SEQ)
    ab = make_antibody_fasta(tmp_path)
    out = tmp_path / "out"

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--sliding",
        "-k", str(k),
        "--model", "fast",
    ])
    assert res.returncode == 0, res.stderr

    windows = generate_sliding_windows(PEPTIDE_SEQ, k)
    assert len(windows) == len(PEPTIDE_SEQ) - k + 1  # sanity

    files = sorted(out.glob("*.json"))
    assert len(files) == len(windows)

    for i, window_seq in enumerate(windows):
        # 2-digit suffix in filename and name
        fname = out / f"pep1_mab1_{k}mer_window_{i:02d}.json"
        assert fname.exists(), f"missing {fname}"
        spec = json.loads(fname.read_text())
        assert spec["name"] == f"pep1_mab1_{k}mer_window_{i:02d}"
        assert spec["sequences"][0]["sequence"] == window_seq
        assert spec["sequences"][1]["sequence"] == HEAVY_SEQ
        assert spec["sequences"][2]["sequence"] == LIGHT_SEQ


def test_sliding_shorter_than_k(tmp_path):
    k = 30  # longer than the 20aa peptide
    pep = make_peptide_fasta(tmp_path, name="pep1", seq=PEPTIDE_SEQ)
    ab = make_antibody_fasta(tmp_path)
    out = tmp_path / "out"

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--sliding",
        "-k", str(k),
        "--model", "fast",
    ])
    assert res.returncode == 0, res.stderr

    files = sorted(out.glob("*.json"))
    assert len(files) == 1
    spec = json.loads(files[0].read_text())
    assert files[0].name == f"pep1_mab1_{k}mer_window_00.json"
    assert spec["sequences"][0]["sequence"] == PEPTIDE_SEQ


def _build_msa_fixture(tmp_path, chain_dirs):
    """Create an msa index + staged msa dir.

    chain_dirs: list of (lookup_key, dirname, create_a3m_bool)
    Returns (index_path, staged_dir_path).
    """
    staged = tmp_path / "staged"
    staged.mkdir()
    index = {}
    for key, dirname, create in chain_dirs:
        d = staged / dirname
        d.mkdir()
        if create:
            (d / "unpaired_msa.a3m").write_text(">q\nACDE\n")
            (d / "paired_msa.a3m").write_text(">q\nACDE\n")
        # msa_dir in the index points at a "published" location whose
        # .name matches the staged subdir name (mirrors real index).
        index[key] = {
            "msa_dir": f"/published/{dirname}",
            "unpaired_msa": f"/published/{dirname}/unpaired_msa.a3m",
            "paired_msa": f"/published/{dirname}/paired_msa.a3m",
            "templates_dir": f"/published/{dirname}/templates",
            "has_templates": False,
            "num_templates": 0,
        }
    index_path = tmp_path / "msa_index.json"
    index_path.write_text(json.dumps(index))
    return str(index_path), str(staged)


def test_full_model_attaches_msas(tmp_path):
    pep = make_peptide_fasta(tmp_path, name="pep1")
    ab = make_antibody_fasta(tmp_path, heavy_id="mab1_heavy", light_id="mab1_light")
    out = tmp_path / "out"

    # peptide + heavy have MSAs; light does NOT (no a3m staged)
    index_path, staged = _build_msa_fixture(tmp_path, [
        ("pep1", "pep1_dir", True),
        ("mab1_heavy", "heavy_dir", True),
        ("mab1_light", "light_dir", False),
    ])

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--model", "full",
        "--msa-index", index_path,
        "--msa-dir", staged,
    ])
    assert res.returncode == 0, res.stderr

    spec = json.loads((out / "pep1_mab1.json").read_text())
    assert "msas" in spec
    msa_chains = {m["chain_id"]: m for m in spec["msas"]}
    # A (peptide) and B (heavy) present, C (light) absent
    assert set(msa_chains) == {"A", "B"}
    for m in spec["msas"]:
        assert "a3m_path" in m
        assert os.path.isabs(m["a3m_path"])
        assert "sequences" not in m  # exactly one source
    # Uses the published unpaired_msa path from the index
    assert msa_chains["A"]["a3m_path"] == "/published/pep1_dir/unpaired_msa.a3m"
    assert msa_chains["B"]["a3m_path"] == "/published/heavy_dir/unpaired_msa.a3m"


def test_fast_model_ignores_msa_index(tmp_path):
    pep = make_peptide_fasta(tmp_path, name="pep1")
    ab = make_antibody_fasta(tmp_path)
    out = tmp_path / "out"
    index_path, staged = _build_msa_fixture(tmp_path, [
        ("pep1", "pep1_dir", True),
        ("mab1_heavy", "heavy_dir", True),
        ("mab1_light", "light_dir", True),
    ])

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--model", "fast",
        "--msa-index", index_path,
        "--msa-dir", staged,
    ])
    assert res.returncode == 0, res.stderr
    spec = json.loads((out / "pep1_mab1.json").read_text())
    assert "msas" not in spec, "fast model must never emit msas even with index"


@pytest.mark.parametrize("heavy_id,light_id", [
    ("mab1_heavy", "mab1_light"),
    ("mab1_HC", "mab1_LC"),
])
def test_chain_detection_variants(tmp_path, heavy_id, light_id):
    pep = make_peptide_fasta(tmp_path, name="pep1")
    ab = make_antibody_fasta(tmp_path, heavy_id=heavy_id, light_id=light_id)
    out = tmp_path / "out"

    res = run_script([
        "--peptide-fasta", pep,
        "--antibody-fasta", ab,
        "--antibody-name", "mab1",
        "--output-dir", str(out),
        "--model", "fast",
    ])
    assert res.returncode == 0, res.stderr
    spec = json.loads((out / "pep1_mab1.json").read_text())
    assert spec["sequences"][1]["sequence"] == HEAVY_SEQ
    assert spec["sequences"][2]["sequence"] == LIGHT_SEQ
