#!/usr/bin/env python3
"""Fold an ESMFold2 JobSpec across multiple seeds, AF3-compatible outputs.

Loads the ESMFold2 model once, folds N diffusion samples in batched GPU forward
passes (via num_diffusion_samples), keeps the best sample by the
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
        # Intentionally chain_iptm[0] (peptide interface), NOT AF3's global
        # ranking_score: PAIRIS ranks on chain_iptm[0]. Do not "fix" this.
        "ranking_score": chain_iptm[0],
        "best_seed": best_seed,
    }


def ensure_atom_site_occupancy(mmcif: str) -> str:
    """Add an ``_atom_site.occupancy`` column (constant 1.00) if absent.

    The Rosetta stage reads structures with Biopython's ``MMCIFParser``, which
    requires ``_atom_site.occupancy``. AF3's mmCIF includes it, but ESMFold2's
    ``to_mmcif()`` omits it (it writes ``B_iso_or_equiv`` only), so its CIFs
    fail to parse downstream. Inject a constant-1.00 occupancy column so the
    output stays drop-in compatible with the AF3 path. Idempotent.
    """
    lines = mmcif.splitlines()
    header_idx = [i for i, ln in enumerate(lines)
                  if ln.strip().startswith("_atom_site.")]
    if not header_idx:
        return mmcif  # no atom_site loop to fix
    if any(lines[i].strip() == "_atom_site.occupancy" for i in header_idx):
        return mmcif  # already present

    # Append occupancy as the final loop column...
    out = list(lines)
    out.insert(header_idx[-1] + 1, "_atom_site.occupancy")
    # ...and a matching value at the end of every atom record.
    result = []
    for ln in out:
        first = ln.split(maxsplit=1)[0] if ln.split() else ""
        result.append(ln.rstrip() + " 1.00" if first in ("ATOM", "HETATM")
                      else ln)
    text = "\n".join(result)
    return text + "\n" if mmcif.endswith("\n") else text


def _rank_key(record) -> tuple:
    """Seed-ranking key: peptide interface chain_iptm[0], tie-break by ptm.

    Single source of truth for the ranking rule, shared by ``pick_best`` and
    the running-best comparison in ``run`` so the two can never diverge.
    """
    return (record["chain_iptm"][0], record["ptm"])


def pick_best(seed_records: list) -> dict:
    """Pick the best seed record: max chain_iptm[0], tie-break by ptm."""
    return max(seed_records, key=_rank_key)


def plan_sample_batches(num_samples: int, per_pass: int) -> list:
    """Split ``num_samples`` diffusion samples into per-forward-pass chunk sizes.

    ``per_pass <= 0`` (or ``>= num_samples``) folds all samples in a single pass;
    otherwise full chunks of ``per_pass`` plus a final remainder. Lets a caller
    cap how many diffusion samples share one GPU forward pass (a VRAM valve)
    while still producing ``num_samples`` total. Pure / GPU-free (unit-tested).
    """
    if num_samples <= 0:
        return []
    if not per_pass or per_pass <= 0 or per_pass >= num_samples:
        return [num_samples]
    full, rem = divmod(num_samples, per_pass)
    chunks = [per_pass] * full
    if rem:
        chunks.append(rem)
    return chunks


# ---------------------------------------------------------------------------
# GPU path (imports torch/esm/transformers lazily; not exercised in tests)
# ---------------------------------------------------------------------------

def load_model(model: str):
    """Load the ESMFold2 model once and move it to GPU in eval mode."""
    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for ESMFold2 inference")
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


def fold_samples(model, spi, num_loops: int, num_sampling_steps: int,
                 num_samples: int, samples_per_pass: int = 0) -> list:
    """Fold ``num_samples`` diffusion samples of one input, batched on the GPU.

    Instead of ``num_samples`` separate single-sample ``fold()`` calls (each
    ``num_diffusion_samples=1``, i.e. batch-of-one, which barely occupies the
    GPU for small peptides), request ``num_diffusion_samples=k`` per forward
    pass so the model works a full batch at once. ``samples_per_pass`` caps
    ``k`` for VRAM (``<=0`` folds all samples in one pass). Returns a flat list
    of per-sample result objects (the same objects the old per-seed loop
    produced), so downstream ranking / output logic is unchanged.
    """
    import torch
    from esm.models.esmfold2 import ESMFold2InputBuilder

    builder = ESMFold2InputBuilder()
    results = []
    for chunk_idx, k in enumerate(plan_sample_batches(num_samples,
                                                      samples_per_pass)):
        with torch.inference_mode():
            out = builder.fold(
                model,
                spi,
                num_loops=num_loops,
                num_sampling_steps=num_sampling_steps,
                num_diffusion_samples=k,
                seed=chunk_idx,
            )
        results.extend(out if isinstance(out, list) else [out])
    return results


def run(spec: dict, base_dir: Path, output_dir, model: str,
        num_seeds: int, samples_per_pass: int = 0) -> int:
    """Load model, fold N samples (batched), write AF3 outputs for the best."""
    outer = spec["name"]

    # Peptide-interface ranking (chain_iptm[0]) and the AF3-compatible outputs
    # assume the peptide is chain A / sequences[0] (the generator always emits it
    # first, and collation reads id == 'A' as the peptide). Fail loudly if a
    # caller ever reorders chains, rather than silently ranking on the wrong one.
    first_id = spec["sequences"][0]["id"]
    first_id = first_id[0] if isinstance(first_id, (list, tuple)) else first_id
    if first_id != "A":
        raise ValueError(
            f"expected peptide as chain 'A' (sequences[0]); got {first_id!r}"
        )

    cif_path, summary_path, all_seeds_path = output_paths(output_dir, outer)
    cif_path.parent.mkdir(parents=True, exist_ok=True)

    config = spec.get("config", {})
    num_loops = config.get("num_loops", 3)
    num_sampling_steps = config.get("num_sampling_steps", 50)

    loaded = load_model(model)
    spi = build_input(spec, base_dir)

    # Fold all `num_seeds` samples in batched GPU forward passes instead of one
    # single-sample fold() per seed. num_diffusion_samples fills the batch
    # dimension, which is what keeps the GPU busy for these small peptide inputs;
    # samples_per_pass caps how many samples share one pass (VRAM valve).
    batches = plan_sample_batches(num_seeds, samples_per_pass)
    t0 = time.perf_counter()
    print(f"Folding {outer!r} ({len(spec['sequences'])} chains): "
          f"{num_seeds} samples in passes of {batches}...")
    results = fold_samples(loaded, spi, num_loops, num_sampling_steps,
                           num_seeds, samples_per_pass)
    print(f"  {len(results)} samples folded in {time.perf_counter() - t0:.1f}s")

    seed_records = []
    best_record = None
    best_mmcif = None
    for i, result in enumerate(results):
        pair_chains_iptm = result.pair_chains_iptm.detach().cpu().tolist()
        chain_iptm = compute_chain_iptm(pair_chains_iptm)
        record = {
            "seed": i,
            "chain_iptm": chain_iptm,
            "ptm": float(result.ptm),
            "iptm": float(result.iptm),
        }
        seed_records.append(record)

        # Track the running best by the shared ranking rule; capture its mmCIF
        # now so non-best results can be dropped.
        if best_record is None or _rank_key(record) > _rank_key(best_record):
            best_record = record
            best_mmcif = result.complex.to_mmcif()

    cif_path.write_text(ensure_atom_site_occupancy(best_mmcif))
    print(f"Wrote {cif_path}")

    best = best_record
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
                   help="Number of diffusion samples to fold, keeping the best "
                        "(default: 5); batched into GPU forward passes.")
    p.add_argument("--samples-per-pass", type=int, default=0,
                   help="Max diffusion samples per GPU forward pass "
                        "(0 = all --num-seeds at once). Lower to cap VRAM; "
                        "higher = better GPU utilization. Set to 1 to reproduce "
                        "the legacy one-sample-at-a-time behaviour.")
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
        samples_per_pass=args.samples_per_pass,
    )


if __name__ == "__main__":
    sys.exit(main())
