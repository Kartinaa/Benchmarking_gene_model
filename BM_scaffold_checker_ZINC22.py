#!/usr/bin/env python3
import os, sys, argparse, gzip, pickle, time
from typing import Set, Tuple
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# ---------- Optional fast compression (uses zstandard if available) ----------
try:
    import zstandard as zstd
    HAVE_ZSTD = True
except Exception:
    HAVE_ZSTD = False

def dump_set(path: str, s: Set[str]) -> None:
    if HAVE_ZSTD and path.endswith(".pkl.zst"):
        cctx = zstd.ZstdCompressor(level=10)
        with open(path, "wb") as f:
            with cctx.stream_writer(f) as zf:
                pickle.dump(s, zf, protocol=5)
    elif path.endswith(".pkl.gz"):
        with gzip.open(path, "wb", compresslevel=6) as f:
            pickle.dump(s, f, protocol=5)
    else:
        with open(path, "wb") as f:
            pickle.dump(s, f, protocol=5)

def load_set(path: str) -> Set[str]:
    if HAVE_ZSTD and path.endswith(".pkl.zst"):
        dctx = zstd.ZstdDecompressor()
        with open(path, "rb") as f:
            with dctx.stream_reader(f) as zf:
                return pickle.load(zf)
    elif path.endswith(".pkl.gz"):
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    else:
        with open(path, "rb") as f:
            return pickle.load(f)

def build_scaffold_set(smi_path: str, n_hint: int = 0) -> Set[str]:
    # Reads first token (SMILES) per line; ignores blanks
    seen = set() if n_hint <= 0 else set()
    with open(smi_path, "r") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line: 
                continue
            smiles = line.split()[0]
            seen.add(smiles)
            if ln % 5_000_000 == 0:
                print(f"[load] processed {ln:,} lines...", file=sys.stderr)
    return seen

def smiles_to_bmscaffold(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    scaf = MurckoScaffold.GetScaffoldForMol(mol)
    if not scaf or scaf.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaf, isomericSmiles=True)

def check_file(path: str, scaffold_set: Set[str]) -> Tuple[float,int,int]:
    total = 0
    matched = 0
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            smi = parts[0]
            bms = smiles_to_bmscaffold(smi)
            if bms:
                total += 1
                if bms in scaffold_set:
                    matched += 1
    pct = (matched / total * 100) if total else 0.0
    return pct, matched, total

def main():
    ap = argparse.ArgumentParser(description="Load scaffold .smi once, then interactively check many .smi files.")
    ap.add_argument("scaffold_smi", help="Huge scaffold .smi (53 GB)")
    ap.add_argument("--cache", default=None, help="Path to cache (e.g., scaffolds.pkl.zst or .pkl.gz). If exists, will load; else will build and save.")
    ap.add_argument("--no-cache-save", action="store_true", help="Do not save cache after building.")
    ap.add_argument("--once", metavar="FILE", help="Check a single .smi and exit (no interactive loop).")
    args = ap.parse_args()

    scaffold_set: Set[str] | None = None

    # Try cache first
    if args.cache and os.path.exists(args.cache):
        t0 = time.time()
        print(f"[cache] loading {args.cache} ...", file=sys.stderr)
        scaffold_set = load_set(args.cache)
        print(f"[cache] loaded {len(scaffold_set):,} entries in {time.time()-t0:.1f}s", file=sys.stderr)

    # Build if needed
    if scaffold_set is None:
        t0 = time.time()
        print(f"[build] reading {args.scaffold_smi} ...", file=sys.stderr)
        raw_set = build_scaffold_set(args.scaffold_smi)
        print(f"[build] {len(raw_set):,} uniques (raw SMILES) in {time.time()-t0:.1f}s", file=sys.stderr)
        # Note: your source file is *already* BM scaffolds; if so, you can keep it as-is.
        # If you want to canonicalize further, do it here.
        scaffold_set = raw_set
        if args.cache and not args.no_cache_save:
            t1 = time.time()
            print(f"[cache] saving to {args.cache} ...", file=sys.stderr)
            dump_set(args.cache, scaffold_set)
            print(f"[cache] saved in {time.time()-t1:.1f}s", file=sys.stderr)

    # Single-shot mode
    if args.once:
        pct, matched, total = check_file(args.once, scaffold_set)
        print(f"{args.once}\t{pct:.2f}%\tmatched={matched}\ttotal={total}")
        return

    # Interactive loop
    print("\nEnter .smi filepath(s) to evaluate (one per line). Type 'quit' to exit.")
    while True:
        try:
            path = input("file> ").strip()
        except EOFError:
            break
        if not path:
            continue
        if path.lower() in {"q", "quit", "exit"}:
            break
        if not os.path.exists(path):
            print(f"[warn] not found: {path}")
            continue
        t0 = time.time()
        pct, matched, total = check_file(path, scaffold_set)
        print(f"{path}\t{pct:.2f}%\tmatched={matched}\ttotal={total}\t({time.time()-t0:.1f}s)")
    print("bye!")

if __name__ == "__main__":
    main()

