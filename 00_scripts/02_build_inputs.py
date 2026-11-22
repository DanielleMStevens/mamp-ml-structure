"""
Build Boltz (1 and 2) and RF3 inputs: one receptor-epitope complex per Model_Building row.

python 00_scripts/02_build_inputs.py \
    SERK_coreceptor.tsv \
    03_out_data/lrr_domain_sequences.fasta \
    04_receptor_only_msas \
    06_input_files_fold_primary_only

Arguments may be omitted (defaults above). Relative arguments resolve from the
repo root, so the script can be launched from any directory. MSA paths written
into the Boltz YAMLs and rf3.json are relative to the repo root, so run Boltz /
RF3 from the repo root on whichever machine holds the repo.
"""
import csv, json, sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULTS = ["SERK_coreceptor.tsv", "03_out_data/lrr_domain_sequences.fasta",
            "04_receptor_only_msas", "06_input_files_fold_primary_only"]


def from_root(p):
    p = Path(p).expanduser()
    return p if p.is_absolute() else (REPO_ROOT / p).resolve()


def portable(p):  # repo-relative when inside the repo, absolute otherwise
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(p)


args = sys.argv[1:5] + DEFAULTS[len(sys.argv[1:5]):]
tsv, fasta, msa_dir, out = map(from_root, args)
(out / "boltz").mkdir(parents=True, exist_ok=True)


def safe(s):  # ColabFold's header -> filename rule
    return "".join(c if c.isalnum() or c in "_.-" else "_" for c in s)


# receptor locus ID (2nd '|' field of each header) -> (header, LRR domain sequence)
domains = {}
for block in open(fasta).read().split(">")[1:]:
    header, *seq = block.strip().split("\n")
    domains.setdefault(safe(header.split("|")[1]), []).append((header, "".join(seq)))

rf3, manifest = [], [["pair_id", "row_id", "Receptor_locus_ID", "Ligand", "a3m", "status"]]
for row in csv.DictReader(open(tsv), delimiter="\t"):
    if row["Source_sheet"] != "Model_Building":
        continue
    locus = row["Receptor_locus_ID"]
    pair_id = f"r{int(row['row_id']):05d}_{safe(locus)}_{safe(row['Ligand'])}"
    hits = domains.get(safe(locus), [])

    if len(hits) != 1:
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], "",
                         "SKIP: not in FASTA" if not hits else "SKIP: duplicate ID in FASTA"])
        continue
    header, receptor = hits[0]
    a3m = msa_dir / f"{safe(header)}.a3m"
    if not a3m.exists():
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], "", "SKIP: no a3m"])
        continue
    if a3m.read_text().split("\n")[1].strip() != receptor:
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], portable(a3m), "SKIP: a3m query != domain"])
        continue

    a3m_rel = portable(a3m)
    epitope = row["Ligand_Sequence"].strip().upper()
    (out / "boltz" / f"{pair_id}.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {receptor}\n      msa: {a3m_rel}\n"
        f"  - protein:\n      id: B\n      sequence: {epitope}\n      msa: empty\n")
    rf3.append({"name": pair_id, "components": [
        {"seq": receptor, "chain_id": "A", "msa_path": a3m_rel},
        {"seq": epitope, "chain_id": "B"}]})
    manifest.append([pair_id, row["row_id"], locus, row["Ligand"], a3m_rel, "built"])

(out / "rf3.json").write_text(json.dumps(rf3, indent=1))
csv.writer(open(out / "manifest.tsv", "w"), delimiter="\t").writerows(manifest)

skipped = [m for m in manifest[1:] if m[5] != "built"]
print(f"built {len(rf3)}, skipped {len(skipped)} -> {portable(out)}")
for m in skipped:
    print("  ", m[2], m[5])
