"""
Build Boltz, RF3, ColabFold/AF2, and AlphaFold3 inputs per Model_Building row.

Modes (mutually exclusive):
  --primary_only  receptor (A) + peptide (B)
                  -> 06_input_files_fold_primary_only/
  --complex       receptor (A) + peptide (B) + co-receptor/SERK (C)
                  -> 07_input_files_complex/

Existing outputs are skipped independently per folding mode (default):
  boltz       boltz/{pair_id}.yaml
  rf3         entry in rf3.json (file rewritten only if entries are added)
  alphafold   colabfold/{pair_id}.fasta + colabfold_msa/{pair_id}.a3m
              (if either AF file is missing, both are rewritten)
  alphafold3  alphafold3/{pair_id}.json

manifest.tsv lists only modes skipped because they already exist
(header-only when nothing was skipped).

Examples:
  python 00_scripts/02_build_inputs.py --primary_only
  python 00_scripts/02_build_inputs.py --complex

Relative paths resolve from the repo root. MSA paths written into Boltz/RF3
inputs are repo-relative; run those tools from the repo root. AlphaFold3
unpairedMsaPath values are relative to each JSON file (AF3 dialect v2).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_TSV = "SERK_coreceptor.tsv"
DEFAULT_RECEPTOR_FASTA = "03_out_data/lrr_domain_sequences.fasta"
DEFAULT_RECEPTOR_MSA = "04_receptor_only_msas"
DEFAULT_SERK_FASTA = "03_out_data/serk_ectodomain_sequences.fasta"
DEFAULT_SERK_MSA = "05_coreceptor_msas"
DEFAULT_PRIMARY_OUT = "06_input_files_fold_primary_only"
DEFAULT_COMPLEX_OUT = "07_input_files_complex"

MANIFEST_STATUS_EXISTS = "skipped: exists"


def from_root(p: str | Path) -> Path:
    p = Path(p).expanduser()
    return p if p.is_absolute() else (REPO_ROOT / p).resolve()


def portable(p: Path) -> str:
    """Repo-relative when inside the repo, absolute otherwise."""
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(p)


def safe(s: str) -> str:
    """ColabFold's header -> filename rule."""
    return "".join(c if c.isalnum() or c in "_.-" else "_" for c in s)


def read_fasta_blocks(path: Path) -> list[tuple[str, str]]:
    blocks = []
    for block in path.read_text().split(">")[1:]:
        header, *seq = block.strip().split("\n")
        blocks.append((header, "".join(seq)))
    return blocks


def a3m_query_sequence(a3m_text: str) -> str:
    return next(
        (ln for ln in a3m_text.splitlines()
         if ln and not ln.startswith(">") and not ln.startswith("#")),
        "",
    )


def strip_a3m_meta(a3m_text: str) -> str:
    lines = a3m_text.replace("\x00", "").splitlines()
    if lines and lines[0].startswith("#") and "\t" in lines[0][1:]:
        lines = lines[1:]
    return "\n".join(lines).strip() + "\n"


def pair_sequences(a3m_chains: list[str], query_sequences: list[str]) -> str:
    """ColabFold pair_sequences for cardinality-1 chains (query row only)."""
    n = len(a3m_chains[0].splitlines())
    paired = [""] * n
    for n_chain, _seq in enumerate(query_sequences):
        lines = a3m_chains[n_chain].splitlines()
        for i, line in enumerate(lines):
            if line.startswith(">"):
                paired[i] += line.replace(">", "\t", 1) if n_chain else line
            else:
                paired[i] += line
    return "\n".join(paired)


def pad_sequences(a3m_chains: list[str], query_sequences: list[str]) -> str:
    """ColabFold pad_sequences: unpaired MSA per chain, gaps for the others."""
    blanks = ["-" * len(seq) for seq in query_sequences]
    out, pos = [], 0
    for n_chain, _seq in enumerate(query_sequences):
        for line in a3m_chains[n_chain].splitlines():
            if not line:
                continue
            if line.startswith(">"):
                out.append(line)
            else:
                out.append("".join(blanks[:pos] + [line] + blanks[pos + 1:]))
        pos += 1
    return "\n".join(out)


def unpaired_multimer_a3m(queries: list[str], msa_texts: list[str | None]) -> str:
    """
    Annotated ColabFold a3m for N chains.

    msa_texts[i] is raw monomer a3m text, or None for single-sequence-only
    (used for the peptide).
    """
    assert len(queries) == len(msa_texts) and len(queries) >= 2
    unpaired = []
    paired = []
    for i, (seq, msa) in enumerate(zip(queries, msa_texts), start=1):
        card = 100 + i
        paired.append(f">{card}\n{seq}\n")
        if msa is None:
            unpaired.append(f">{card}\n{seq}\n")
        else:
            unpaired.append(strip_a3m_meta(msa))
    lengths = ",".join(str(len(s)) for s in queries)
    cards = ",".join("1" for _ in queries)
    header = f"#{lengths}\t{cards}\n"
    return (
        header
        + pair_sequences(paired, queries)
        + "\n"
        + pad_sequences(unpaired, queries)
        + "\n"
    )


def load_cached_a3m(path: Path, cache: dict[Path, str]) -> str:
    if path not in cache:
        cache[path] = path.read_text()
    return cache[path]


def load_rf3_json(path: Path) -> list[dict]:
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        raise SystemExit(f"Expected a JSON list in {path}")
    return data


def write_tsv(path: Path, rows: list[list]) -> None:
    with open(path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(rows)


def relpath_to(path: Path, start: Path) -> str:
    """Path relative to start, suitable for AF3 unpairedMsaPath."""
    return os.path.relpath(path.resolve(), start.resolve()).replace("\\", "/")


def af3_protein(
    chain_id: str,
    sequence: str,
    *,
    msa_path: Path | None,
    json_dir: Path,
) -> dict:
    """
    One AF3 protein chain. Heteromers use separate entries with distinct ids.

    Custom MSA (AF3 docs): unpairedMsaPath + pairedMsa \"\" + templates [].
    No MSA (peptide): unpairedMsa \"\" + pairedMsa \"\" + templates [].
    """
    protein: dict = {"id": chain_id, "sequence": sequence}
    if msa_path is None:
        protein["unpairedMsa"] = ""
        protein["pairedMsa"] = ""
        protein["templates"] = []
    else:
        protein["unpairedMsaPath"] = relpath_to(msa_path, json_dir)
        protein["pairedMsa"] = ""
        protein["templates"] = []
    return {"protein": protein}


def af3_job(name: str, sequences: list[dict]) -> dict:
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 2,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build Boltz / RF3 / ColabFold / AlphaFold3 inputs for "
                    "primary or receptor–peptide–co-receptor complexes."
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--primary_only",
        action="store_true",
        help=f"Receptor + peptide -> {DEFAULT_PRIMARY_OUT}/",
    )
    mode.add_argument(
        "--complex",
        action="store_true",
        help=f"Receptor + peptide + co-receptor -> {DEFAULT_COMPLEX_OUT}/",
    )
    p.add_argument("--tsv", default=DEFAULT_TSV,
                   help=f"Pairing table (default: {DEFAULT_TSV})")
    p.add_argument("--receptor-fasta", default=DEFAULT_RECEPTOR_FASTA,
                   help="Primary receptor LRR-domain FASTA")
    p.add_argument("--receptor-msa-dir", default=DEFAULT_RECEPTOR_MSA,
                   help="Primary receptor monomer a3m directory")
    p.add_argument("--serk-fasta", default=DEFAULT_SERK_FASTA,
                   help="Co-receptor ectodomain FASTA (complex mode)")
    p.add_argument("--serk-msa-dir", default=DEFAULT_SERK_MSA,
                   help="Co-receptor monomer a3m directory (complex mode)")
    p.add_argument(
        "--out",
        default=None,
        help="Output directory (defaults depend on --primary_only / --complex)",
    )
    return p.parse_args(argv)


def build_primary(args: argparse.Namespace) -> None:
    tsv = from_root(args.tsv)
    fasta = from_root(args.receptor_fasta)
    msa_dir = from_root(args.receptor_msa_dir)
    out = from_root(args.out or DEFAULT_PRIMARY_OUT)

    boltz_dir = out / "boltz"
    cf_fa_dir = out / "colabfold"
    cf_msa_dir = out / "colabfold_msa"
    af3_dir = out / "alphafold3"
    rf3_path = out / "rf3.json"
    for d in (boltz_dir, cf_fa_dir, cf_msa_dir, af3_dir):
        d.mkdir(parents=True, exist_ok=True)

    domains: dict[str, list[tuple[str, str]]] = {}
    for header, seq in read_fasta_blocks(fasta):
        parts = header.split("|")
        if len(parts) < 2:
            continue
        domains.setdefault(safe(parts[1]), []).append((header, seq))

    rf3 = load_rf3_json(rf3_path)
    rf3_index = {entry["name"]: i for i, entry in enumerate(rf3)}
    rf3_dirty = False

    # manifest.tsv: only modes skipped because outputs already exist
    manifest = [["pair_id", "row_id", "Receptor_locus_ID", "Ligand",
                 "mode", "status"]]
    a3m_cache: dict[Path, str] = {}
    n_boltz = n_rf3 = n_af = n_af3 = 0
    errors: list[str] = []

    for row in csv.DictReader(open(tsv), delimiter="\t"):
        if row["Source_sheet"] != "Model_Building":
            continue
        locus = row["Receptor_locus_ID"]
        ligand = row["Ligand"]
        pair_id = f"r{int(row['row_id']):05d}_{safe(locus)}_{safe(ligand)}"
        hits = domains.get(safe(locus), [])

        if len(hits) != 1:
            reason = ("not in FASTA" if not hits else "duplicate ID in FASTA")
            errors.append(f"{pair_id}: {reason}")
            continue

        header, receptor = hits[0]
        a3m = msa_dir / f"{safe(header)}.a3m"
        if not a3m.exists():
            errors.append(f"{pair_id}: no receptor a3m")
            continue

        a3m_text = load_cached_a3m(a3m, a3m_cache)
        if a3m_query_sequence(a3m_text) != receptor:
            errors.append(f"{pair_id}: a3m query != domain")
            continue

        a3m_rel = portable(a3m)
        epitope = row["Ligand_Sequence"].strip().upper()
        yaml_path = boltz_dir / f"{pair_id}.yaml"
        fa_path = cf_fa_dir / f"{pair_id}.fasta"
        cf_a3m_path = cf_msa_dir / f"{pair_id}.a3m"
        af3_path = af3_dir / f"{pair_id}.json"
        base = [pair_id, row["row_id"], locus, ligand]

        # --- boltz ---
        if yaml_path.exists():
            manifest.append(base + ["boltz", MANIFEST_STATUS_EXISTS])
        else:
            yaml_path.write_text(
                "version: 1\nsequences:\n"
                f"  - protein:\n      id: A\n      sequence: {receptor}\n      msa: {a3m_rel}\n"
                f"  - protein:\n      id: B\n      sequence: {epitope}\n      msa: empty\n"
            )
            n_boltz += 1

        # --- rf3 (only add missing pair entries) ---
        rf3_entry = {
            "name": pair_id,
            "components": [
                {"seq": receptor, "chain_id": "A", "msa_path": a3m_rel},
                {"seq": epitope, "chain_id": "B"},
            ],
        }
        if pair_id in rf3_index:
            manifest.append(base + ["rf3", MANIFEST_STATUS_EXISTS])
        else:
            rf3.append(rf3_entry)
            rf3_index[pair_id] = len(rf3) - 1
            rf3_dirty = True
            n_rf3 += 1

        # --- alphafold (fasta + a3m as a unit) ---
        if fa_path.exists() and cf_a3m_path.exists():
            manifest.append(base + ["alphafold", MANIFEST_STATUS_EXISTS])
        else:
            fa_path.write_text(f">{pair_id}\n{receptor}:{epitope}\n")
            cf_a3m_path.write_text(
                unpaired_multimer_a3m([receptor, epitope], [a3m_text, None])
            )
            n_af += 1

        # --- alphafold3 ---
        if af3_path.exists():
            manifest.append(base + ["alphafold3", MANIFEST_STATUS_EXISTS])
        else:
            af3_path.write_text(json.dumps(af3_job(pair_id, [
                af3_protein("A", receptor, msa_path=a3m, json_dir=af3_dir),
                af3_protein("B", epitope, msa_path=None, json_dir=af3_dir),
            ]), indent=2) + "\n")
            n_af3 += 1

    if rf3_dirty:
        rf3_path.write_text(json.dumps(rf3, indent=1))
    write_tsv(out / "manifest.tsv", manifest)

    n_skipped = len(manifest) - 1
    print(f"primary_only -> {portable(out)}")
    print(f"  wrote   boltz={n_boltz}  rf3={n_rf3}  alphafold={n_af}  "
          f"alphafold3={n_af3}")
    print(f"  skipped {n_skipped} mode(s) already present -> manifest.tsv")
    if rf3_dirty:
        print(f"  updated rf3.json ({len(rf3)} entries)")
    else:
        print("  rf3.json unchanged")
    for msg in errors[:30]:
        print("  ERROR", msg)
    if len(errors) > 30:
        print(f"  ... and {len(errors) - 30} more errors")


def build_complex(args: argparse.Namespace) -> None:
    tsv = from_root(args.tsv)
    receptor_fasta = from_root(args.receptor_fasta)
    receptor_msa_dir = from_root(args.receptor_msa_dir)
    serk_fasta = from_root(args.serk_fasta)
    serk_msa_dir = from_root(args.serk_msa_dir)
    out = from_root(args.out or DEFAULT_COMPLEX_OUT)

    boltz_dir = out / "boltz"
    cf_fa_dir = out / "colabfold"
    cf_msa_dir = out / "colabfold_msa"
    af3_dir = out / "alphafold3"
    rf3_path = out / "rf3.json"
    for d in (boltz_dir, cf_fa_dir, cf_msa_dir, af3_dir):
        d.mkdir(parents=True, exist_ok=True)

    domains: dict[str, list[tuple[str, str]]] = {}
    for header, seq in read_fasta_blocks(receptor_fasta):
        parts = header.split("|")
        if len(parts) < 2:
            continue
        domains.setdefault(safe(parts[1]), []).append((header, seq))

    serk_ecto: dict[str, str] = {}
    for header, seq in read_fasta_blocks(serk_fasta):
        serk_id = header.split()[0]
        if serk_id in serk_ecto and serk_ecto[serk_id] != seq:
            raise SystemExit(f"Duplicate conflicting SERK ectodomain for {serk_id}")
        serk_ecto[serk_id] = seq

    rf3 = load_rf3_json(rf3_path)
    rf3_index = {entry["name"]: i for i, entry in enumerate(rf3)}
    rf3_dirty = False

    manifest = [["pair_id", "row_id", "Receptor_locus_ID", "Ligand", "SERK_ID",
                 "mode", "status"]]
    missing_serk = [[
        "pair_id", "row_id", "Receptor_locus_ID", "Ligand", "SERK_ID", "reason",
    ]]
    a3m_cache: dict[Path, str] = {}
    n_boltz = n_rf3 = n_af = n_af3 = 0
    errors: list[str] = []

    for row in csv.DictReader(open(tsv), delimiter="\t"):
        if row["Source_sheet"] != "Model_Building":
            continue
        locus = row["Receptor_locus_ID"]
        ligand = row["Ligand"]
        pair_id = f"r{int(row['row_id']):05d}_{safe(locus)}_{safe(ligand)}"
        serk_id = (row.get("SERK_ID") or "").strip()

        hits = domains.get(safe(locus), [])
        if len(hits) != 1:
            reason = ("not in FASTA" if not hits else "duplicate ID in FASTA")
            errors.append(f"{pair_id}: {reason}")
            continue

        header, receptor = hits[0]
        rec_a3m = receptor_msa_dir / f"{safe(header)}.a3m"
        if not rec_a3m.exists():
            errors.append(f"{pair_id}: no receptor a3m")
            continue

        rec_a3m_text = load_cached_a3m(rec_a3m, a3m_cache)
        if a3m_query_sequence(rec_a3m_text) != receptor:
            errors.append(f"{pair_id}: receptor a3m query != domain")
            continue

        if not serk_id or serk_id.upper() == "NA":
            reason = "SKIP: no SERK_ID"
            missing_serk.append([pair_id, row["row_id"], locus, ligand,
                                 serk_id, reason])
            continue

        coreceptor = serk_ecto.get(serk_id)
        if not coreceptor:
            reason = "SKIP: SERK_ID not in ectodomain FASTA"
            missing_serk.append([pair_id, row["row_id"], locus, ligand,
                                 serk_id, reason])
            continue

        serk_a3m = serk_msa_dir / f"{serk_id}.a3m"
        if not serk_a3m.exists():
            serk_a3m_safe = serk_msa_dir / f"{safe(serk_id)}.a3m"
            if serk_a3m_safe.exists():
                serk_a3m = serk_a3m_safe
            else:
                reason = "SKIP: no SERK a3m"
                missing_serk.append([pair_id, row["row_id"], locus, ligand,
                                     serk_id, reason])
                continue

        serk_a3m_text = load_cached_a3m(serk_a3m, a3m_cache)
        if a3m_query_sequence(serk_a3m_text) != coreceptor:
            reason = "SKIP: SERK a3m query != ectodomain"
            missing_serk.append([pair_id, row["row_id"], locus, ligand,
                                 serk_id, reason])
            continue

        epitope = row["Ligand_Sequence"].strip().upper()
        rec_a3m_rel = portable(rec_a3m)
        serk_a3m_rel = portable(serk_a3m)
        yaml_path = boltz_dir / f"{pair_id}.yaml"
        fa_path = cf_fa_dir / f"{pair_id}.fasta"
        cf_a3m_path = cf_msa_dir / f"{pair_id}.a3m"
        af3_path = af3_dir / f"{pair_id}.json"
        base = [pair_id, row["row_id"], locus, ligand, serk_id]

        # --- boltz ---
        if yaml_path.exists():
            manifest.append(base + ["boltz", MANIFEST_STATUS_EXISTS])
        else:
            yaml_path.write_text(
                "version: 1\nsequences:\n"
                f"  - protein:\n      id: A\n      sequence: {receptor}\n      msa: {rec_a3m_rel}\n"
                f"  - protein:\n      id: B\n      sequence: {epitope}\n      msa: empty\n"
                f"  - protein:\n      id: C\n      sequence: {coreceptor}\n      msa: {serk_a3m_rel}\n"
            )
            n_boltz += 1

        # --- rf3 ---
        rf3_entry = {
            "name": pair_id,
            "components": [
                {"seq": receptor, "chain_id": "A", "msa_path": rec_a3m_rel},
                {"seq": epitope, "chain_id": "B"},
                {"seq": coreceptor, "chain_id": "C", "msa_path": serk_a3m_rel},
            ],
        }
        if pair_id in rf3_index:
            manifest.append(base + ["rf3", MANIFEST_STATUS_EXISTS])
        else:
            rf3.append(rf3_entry)
            rf3_index[pair_id] = len(rf3) - 1
            rf3_dirty = True
            n_rf3 += 1

        # --- alphafold (fasta + a3m as a unit) ---
        if fa_path.exists() and cf_a3m_path.exists():
            manifest.append(base + ["alphafold", MANIFEST_STATUS_EXISTS])
        else:
            fa_path.write_text(
                f">{pair_id}\n{receptor}:{epitope}:{coreceptor}\n"
            )
            cf_a3m_path.write_text(
                unpaired_multimer_a3m(
                    [receptor, epitope, coreceptor],
                    [rec_a3m_text, None, serk_a3m_text],
                )
            )
            n_af += 1

        # --- alphafold3 ---
        if af3_path.exists():
            manifest.append(base + ["alphafold3", MANIFEST_STATUS_EXISTS])
        else:
            af3_path.write_text(json.dumps(af3_job(pair_id, [
                af3_protein("A", receptor, msa_path=rec_a3m, json_dir=af3_dir),
                af3_protein("B", epitope, msa_path=None, json_dir=af3_dir),
                af3_protein("C", coreceptor, msa_path=serk_a3m, json_dir=af3_dir),
            ]), indent=2) + "\n")
            n_af3 += 1

    if rf3_dirty:
        rf3_path.write_text(json.dumps(rf3, indent=1))
    write_tsv(out / "manifest.tsv", manifest)
    write_tsv(out / "missing_serk_msa.tsv", missing_serk)

    n_skipped = len(manifest) - 1
    n_missing = len(missing_serk) - 1
    print(f"complex -> {portable(out)}")
    print(f"  wrote   boltz={n_boltz}  rf3={n_rf3}  alphafold={n_af}  "
          f"alphafold3={n_af3}")
    print(f"  skipped {n_skipped} mode(s) already present -> manifest.tsv")
    if rf3_dirty:
        print(f"  updated rf3.json ({len(rf3)} entries)")
    else:
        print("  rf3.json unchanged")
    print(f"  missing_serk_msa.tsv {n_missing} rows")
    for msg in errors[:30]:
        print("  ERROR", msg)
    if len(errors) > 30:
        print(f"  ... and {len(errors) - 30} more errors")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.primary_only:
        build_primary(args)
    else:
        build_complex(args)


if __name__ == "__main__":
    main()
