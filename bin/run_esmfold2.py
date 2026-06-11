#!/usr/bin/env python3
"""Fold an ESMFold2 JobSpec across multiple seeds, AF3-compatible outputs.

Loads the ESMFold2 model once, folds seeds 0..N-1, keeps the best seed by the
peptide-interface confidence, and writes outputs in the layout PAIRIS's existing
AlphaFold3 (AF3) downstream stages consume unchanged. For an input JSON whose
``name`` is OUTER, with INNER = OUTER.lower():

    <output_dir>/<INNER>/<INNER>_model.cif
    <output_dir>/<INNER>/<INNER>_summary_confidences.json
    <output_dir>/<INNER>/<INNER>_all_seeds.json   (diagnostics)

The reported peptide-interface metric is the AF3 analog of ``chain_iptm[0]``:
chain A = peptide, so ``chain_iptm[0]`` = mean of the peptide's inter-chain
expected-TM values, derived from the model's ``pair_chains_iptm`` tensor.

Heavy imports (torch / esm / transformers) live inside the functions that need
them so the pure helpers below can be imported and unit-tested without a GPU.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

# Same mapping as the reference runner.py.
MODEL_IDS = {
    "fast": "biohub/ESMFold2-Fast",
    "full": "biohub/ESMFold2",
}


# ---------------------------------------------------------------------------
# Pure, GPU-free helpers (unit-tested without torch/esm/transformers)
# ---------------------------------------------------------------------------

def compute_chain_iptm(pair_chains_iptm) -> list:
    """Per-chain interface means from a square pair-chains iPTM matrix.

    ``chain_iptm[i] = mean over j != i of m[i][j]``. Accepts a nested
    list/array (converts to plain lists internally). For a single chain there
    are no off-diagonal entries, so returns ``[0.0]``.
    """
    # Convert to nested lists so we don't depend on numpy/torch here.
    m = [list(row) for row in pair_chains_iptm]
    n = len(m)
    if n <= 1:
        return [0.0]

    result = []
    for i in range(n):
        off_diag = [float(m[i][j]) for j in range(n) if j != i]
        result.append(sum(off_diag) / len(off_diag))
    return result


def output_paths(output_dir, outer_name: str):
    """Return (cif_path, summary_path, all_seeds_path) for an OUTER name.

    INNER = outer_name.lower(); files nest under <output_dir>/<INNER>/.
    """
    inner = outer_name.lower()
    base = Path(output_dir) / inner
    cif_path = base / f"{inner}_model.cif"
    summary_path = base / f"{inner}_summary_confidences.json"
    all_seeds_path = base / f"{inner}_all_seeds.json"
    return cif_path, summary_path, all_seeds_path


def build_summary_confidences(chain_iptm: list, ptm: float, iptm: float,
                              best_seed: int) -> dict:
    """Build the AF3-compatible summary_confidences dict for the best seed.

    ``chain_iptm`` order is preserved (peptide chain A first). ``ranking_score``
    is the peptide interface value chain_iptm[0], matching how PAIRIS ranks.
    """
    return {
        "chain_iptm": list(chain_iptm),
        "ptm": ptm,
        "iptm": iptm,
        "ranking_score": chain_iptm[0],
        "best_seed": best_seed,
    }


def pick_best(seed_records: list) -> dict:
    """Pick the best seed record: max chain_iptm[0], tie-break by ptm."""
    return max(
        seed_records,
        key=lambda r: (r["chain_iptm"][0], r["ptm"]),
    )


# ---------------------------------------------------------------------------
# GPU path (imports torch/esm/transformers lazily; not exercised in tests)
# ---------------------------------------------------------------------------

def load_model(model: str):
    """Load the ESMFold2 model once and move it to GPU in eval mode."""
    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

    assert torch.cuda.is_available(), "CUDA required for ESMFold2 inference"
    print(f"torch={torch.__version__} cuda={torch.version.cuda}")
    print(f"GPU: {torch.cuda.get_device_name(0)}")

    model_id = MODEL_IDS[model]
    t0 = time.perf_counter()
    print(f"Loading {model_id}...")
    loaded = ESMFold2Model.from_pretrained(model_id).cuda().eval()
    print(f"Model loaded in {time.perf_counter() - t0:.1f}s")
    return loaded


def build_input(spec: dict, base_dir: Path):
    """Build a StructurePredictionInput once from the JobSpec dict.

    Faithful port of esmfold2-test/builder.py for the entity types PAIRIS
    emits (protein chains + optional per-chain MSA). MSA ``a3m_path`` is
    resolved as ``base_dir / a3m_path`` (absolute paths pass through).
    """
    from esm.utils.msa import MSA
    from esm.utils.structure import input_builder as ib

    # Map chain_id -> MSA (only protein chains carry MSAs in PAIRIS).
    msa_by_chain: dict = {}
    for e in spec.get("msas", []):
        chain_id = e["chain_id"]
        a3m_path = e.get("a3m_path")
        if a3m_path is not None:
            msa_by_chain[chain_id] = MSA.from_a3m(
                path=str(base_dir / a3m_path),
                remove_insertions=e.get("remove_insertions", True),
                max_sequences=e.get("max_sequences"),
            )
        else:
            msa_by_chain[chain_id] = MSA.from_sequences(e["sequences"])

    sequences = []
    for s in spec["sequences"]:
        seq_type = s["type"]
        if seq_type != "protein":
            # PAIRIS never emits ligands/RNA/DNA/mods/bonds; be explicit.
            raise ValueError(
                f"run_esmfold2 only supports protein chains; got {seq_type!r}"
            )
        ids = [s["id"]] if isinstance(s["id"], str) else list(s["id"])
        msas_for_entity = [msa_by_chain[cid] for cid in ids if cid in msa_by_chain]
        if len(msas_for_entity) > 1:
            raise ValueError(
                f"multiple MSAs assigned to one multi-chain entity {ids}"
            )
        kwargs: dict = {"id": s["id"], "sequence": s["sequence"]}
        if msas_for_entity:
            kwargs["msa"] = msas_for_entity[0]
        sequences.append(ib.ProteinInput(**kwargs))

    return ib.StructurePredictionInput(sequences=sequences)


def fold_seed(model, spi, num_loops: int, num_sampling_steps: int, seed: int):
    """Fold a single seed and return the model result object."""
    import torch
    from esm.models.esmfold2 import ESMFold2InputBuilder

    with torch.inference_mode():
        return ESMFold2InputBuilder().fold(
            model,
            spi,
            num_loops=num_loops,
            num_sampling_steps=num_sampling_steps,
            seed=seed,
        )


def run(spec: dict, base_dir: Path, output_dir, model: str,
        num_seeds: int) -> int:
    """Load model, fold all seeds, write AF3-compatible outputs for the best."""
    outer = spec["name"]
    cif_path, summary_path, all_seeds_path = output_paths(output_dir, outer)
    cif_path.parent.mkdir(parents=True, exist_ok=True)

    config = spec.get("config", {})
    num_loops = config.get("num_loops", 3)
    num_sampling_steps = config.get("num_sampling_steps", 50)

    loaded = load_model(model)
    spi = build_input(spec, base_dir)

    seed_records = []
    best_record = None
    best_result = None
    for seed in range(num_seeds):
        t0 = time.perf_counter()
        print(f"Folding {outer!r} seed {seed} "
              f"({len(spec['sequences'])} chains)...")
        result = fold_seed(loaded, spi, num_loops, num_sampling_steps, seed)
        print(f"  seed {seed} done in {time.perf_counter() - t0:.1f}s")

        pair_chains_iptm = result.pair_chains_iptm.detach().cpu().tolist()
        chain_iptm = compute_chain_iptm(pair_chains_iptm)
        ptm = float(result.ptm)
        iptm = float(result.iptm)
        record = {
            "seed": seed,
            "chain_iptm": chain_iptm,
            "ptm": ptm,
            "iptm": iptm,
        }
        seed_records.append(record)

        # Track the running best so we can keep its result (cif) without
        # holding every seed's structure in memory.
        if best_record is None or (
            (record["chain_iptm"][0], record["ptm"])
            > (best_record["chain_iptm"][0], best_record["ptm"])
        ):
            best_record = record
            best_result = result

    best = pick_best(seed_records)
    assert best["seed"] == best_record["seed"], "running-best disagreed with pick_best"

    cif_path.write_text(best_result.complex.to_mmcif())
    print(f"Wrote {cif_path}")

    summary = build_summary_confidences(
        chain_iptm=best["chain_iptm"],
        ptm=best["ptm"],
        iptm=best["iptm"],
        best_seed=best["seed"],
    )
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {summary_path}")

    all_seeds = [
        {
            "seed": r["seed"],
            "chain_iptm0": r["chain_iptm"][0],
            "ptm": r["ptm"],
            "iptm": r["iptm"],
        }
        for r in seed_records
    ]
    all_seeds_path.write_text(json.dumps(all_seeds, indent=2))
    print(f"Wrote {all_seeds_path}")

    print(f"Best seed: {best['seed']} (chain_iptm[0]={best['chain_iptm'][0]:.4f})")
    return 0


def main(argv: list = None) -> int:
    p = argparse.ArgumentParser(
        description="Fold an ESMFold2 JobSpec across seeds; write "
                    "AF3-compatible outputs.")
    p.add_argument("--input", required=True, type=Path,
                   help="ESMFold2 JobSpec JSON")
    p.add_argument("--output-dir", type=Path, default=Path("."),
                   help="Output directory (default: .)")
    p.add_argument("--model", choices=["fast", "full"], default="fast",
                   help="ESMFold2 model variant (default: fast)")
    p.add_argument("--num-seeds", type=int, default=5,
                   help="Number of seeds to fold (seeds 0..N-1, default: 5)")
    p.add_argument("--overwrite", action="store_true",
                   help="Replace existing outputs")
    args = p.parse_args(argv)

    if not args.input.exists():
        print(f"ERROR: input not found: {args.input}", file=sys.stderr)
        return 2

    spec = json.loads(args.input.read_text())
    outer = spec["name"]
    cif_path, summary_path, _ = output_paths(args.output_dir, outer)

    # Idempotent resume: skip already-completed work unless --overwrite.
    if not args.overwrite and cif_path.exists() and summary_path.exists():
        print(f"Outputs already exist for {outer!r}; skipping "
              f"(pass --overwrite to replace).")
        return 0

    return run(
        spec=spec,
        base_dir=args.input.parent,
        output_dir=args.output_dir,
        model=args.model,
        num_seeds=args.num_seeds,
    )


if __name__ == "__main__":
    sys.exit(main())
