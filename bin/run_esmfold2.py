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


def input_output_dir(base_output_dir, input_path: Path, batch: bool) -> Path:
    """Per-input output directory.

    In batch mode each input nests under its own stem (``base/<stem>/``) so many
    inputs can share a single ``--output-dir`` without colliding; in single mode
    the base dir is used directly, preserving the legacy one-input layout
    exactly (the Nextflow module passes ``--output-dir <json.baseName>``).
    """
    base = Path(base_output_dir)
    return base / input_path.stem if batch else base


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


def run(loaded, spec: dict, base_dir: Path, output_dir,
        num_seeds: int) -> int:
    """Fold all seeds with a pre-loaded model; write AF3 outputs for the best.

    ``loaded`` is the model returned by :func:`load_model`. It is passed in
    (rather than loaded here) so one process can fold many specs after a single
    model load -- the dominant cost -- instead of reloading it per input.
    """
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

    spi = build_input(spec, base_dir)

    seed_records = []
    best_record = None
    best_mmcif = None
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

        # Track the running best by the shared ranking rule and capture its
        # mmCIF as a host string now, so each seed's GPU result can go out of
        # scope this iteration (peak memory = one live result).
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
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", type=Path,
                     help="A single ESMFold2 JobSpec JSON")
    src.add_argument("--inputs", type=Path, nargs="+",
                     help="Several JobSpec JSONs folded in ONE process: the "
                          "model is loaded once and the specs folded "
                          "sequentially, then each input's outputs nest under "
                          "<output-dir>/<input-stem>/. Amortises the (dominant) "
                          "model-load cost across many small folds.")
    p.add_argument("--output-dir", type=Path, default=Path("."),
                   help="Output directory (default: .)")
    p.add_argument("--model", choices=["fast", "full"], default="fast",
                   help="ESMFold2 model variant (default: fast)")
    p.add_argument("--num-seeds", type=int, default=5,
                   help="Number of seeds to fold (seeds 0..N-1, default: 5)")
    p.add_argument("--overwrite", action="store_true",
                   help="Replace existing outputs")
    args = p.parse_args(argv)

    batch = args.inputs is not None
    inputs = args.inputs if batch else [args.input]

    # Preserve legacy single-input behaviour exactly (including exit code 2).
    if not batch and not args.input.exists():
        print(f"ERROR: input not found: {args.input}", file=sys.stderr)
        return 2

    # Load the model ONCE for the whole (possibly batched) invocation -- this is
    # the expensive step (~tens of seconds), so folding many specs per load is
    # the entire point of --inputs.
    loaded = load_model(args.model)

    n_ok = n_skip = n_fail = 0
    for inp in inputs:
        if not inp.exists():
            print(f"ERROR: input not found: {inp}", file=sys.stderr)
            n_fail += 1
            continue
        out_dir = input_output_dir(args.output_dir, inp, batch)
        try:
            spec = json.loads(inp.read_text())
            outer = spec["name"]
            cif_path, summary_path, _ = output_paths(out_dir, outer)

            # Idempotent resume: skip already-completed work unless --overwrite.
            if not args.overwrite and cif_path.exists() and summary_path.exists():
                print(f"Outputs already exist for {outer!r}; skipping "
                      f"(pass --overwrite to replace).")
                n_skip += 1
                continue

            run(loaded, spec, base_dir=inp.parent, output_dir=out_dir,
                num_seeds=args.num_seeds)
            n_ok += 1
        except Exception as e:  # a single bad spec must not abort a batch
            print(f"ERROR folding {inp}: {e}", file=sys.stderr)
            n_fail += 1

    if batch:
        print(f"Batch complete: {n_ok} folded, {n_skip} skipped, "
              f"{n_fail} failed of {len(inputs)} inputs.")

    # Only hard-fail the task if nothing at all completed; otherwise let the
    # good outputs publish and rely on Nextflow retry + the idempotent skip
    # above to re-attempt just the unfinished specs on the next try.
    return 1 if n_fail and (n_ok + n_skip) == 0 else 0


if __name__ == "__main__":
    sys.exit(main())
