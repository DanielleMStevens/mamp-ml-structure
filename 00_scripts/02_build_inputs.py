"""
Build Boltz, RF3, and ColabFold/AF2 inputs: one receptor-epitope complex per
Model_Building row.

python 00_scripts/02_build_inputs.py \
    SERK_coreceptor.tsv \
    03_out_data/lrr_domain_sequences.fasta \
    04_receptor_only_msas \
    06_input_files_fold_primary_only

Arguments may be omitted (defaults above). Relative arguments resolve from the
repo root. MSA paths in Boltz/RF3 outputs are repo-relative; run those tools
from the repo root.

ColabFold outputs under <out>/:
  colabfold/{pair_id}.fasta     receptor:peptide (MSA server if you use these)
  colabfold_msa/{pair_id}.a3m   unpaired complex MSA reusing the receptor a3m
                                (pass this directory to colabfold_batch; no MSA
                                server). Peptide is single-sequence only.

Example:
  colabfold_batch --model-type alphafold2_multimer_v3 --pair-mode unpaired \\
    06_input_files_fold_primary_only/colabfold_msa \\
    08_output_fold_primary_only/colabfold_af2
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


def safe(s):  # ColabFold's header -> filename rule
    return "".join(c if c.isalnum() or c in "_.-" else "_" for c in s)


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


def complex_unpaired_a3m(receptor: str, epitope: str, receptor_a3m: str) -> str:
    """Annotated ColabFold a3m: receptor MSA unpaired + peptide single sequence."""
    # Drop a leading '#...\t...' meta line if present (monomer a3ms usually lack it).
    lines = receptor_a3m.replace("\x00", "").splitlines()
    if lines and lines[0].startswith("#") and "\t" in lines[0][1:]:
        lines = lines[1:]
    rec_msa = "\n".join(lines).strip() + "\n"
    pep_msa = f">102\n{epitope}\n"
    queries = [receptor, epitope]
    paired = [f">101\n{receptor}\n", f">102\n{epitope}\n"]
    unpaired = [rec_msa, pep_msa]
    header = f"#{len(receptor)},{len(epitope)}\t1,1\n"
    return header + pair_sequences(paired, queries) + "\n" + pad_sequences(unpaired, queries) + "\n"


args = sys.argv[1:5] + DEFAULTS[len(sys.argv[1:5]):]
tsv, fasta, msa_dir, out = map(from_root, args)
boltz_dir = out / "boltz"
cf_fa_dir = out / "colabfold"
cf_msa_dir = out / "colabfold_msa"
for d in (boltz_dir, cf_fa_dir, cf_msa_dir):
    d.mkdir(parents=True, exist_ok=True)


# receptor locus ID (2nd '|' field of each header) -> (header, LRR domain sequence)
domains = {}
for block in open(fasta).read().split(">")[1:]:
    header, *seq = block.strip().split("\n")
    domains.setdefault(safe(header.split("|")[1]), []).append((header, "".join(seq)))

rf3, manifest = [], [["pair_id", "row_id", "Receptor_locus_ID", "Ligand",
                      "a3m", "colabfold_fasta", "colabfold_a3m", "status"]]
# Cache transformed complex a3m text per (receptor a3m path, epitope) is overkill;
# cache per receptor a3m path the raw text instead.
a3m_cache: dict[Path, str] = {}

for row in csv.DictReader(open(tsv), delimiter="\t"):
    if row["Source_sheet"] != "Model_Building":
        continue
    locus = row["Receptor_locus_ID"]
    pair_id = f"r{int(row['row_id']):05d}_{safe(locus)}_{safe(row['Ligand'])}"
    hits = domains.get(safe(locus), [])

    if len(hits) != 1:
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], "", "", "",
                         "SKIP: not in FASTA" if not hits else "SKIP: duplicate ID in FASTA"])
        continue
    header, receptor = hits[0]
    a3m = msa_dir / f"{safe(header)}.a3m"
    if not a3m.exists():
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], "", "", "", "SKIP: no a3m"])
        continue
    if a3m not in a3m_cache:
        a3m_cache[a3m] = a3m.read_text()
    a3m_text = a3m_cache[a3m]
    # Query line is the first sequence after the first '>'
    qseq = next((ln for ln in a3m_text.splitlines()
                 if ln and not ln.startswith(">") and not ln.startswith("#")), "")
    if qseq != receptor:
        manifest.append([pair_id, row["row_id"], locus, row["Ligand"], portable(a3m),
                         "", "", "SKIP: a3m query != domain"])
        continue

    a3m_rel = portable(a3m)
    epitope = row["Ligand_Sequence"].strip().upper()

    (boltz_dir / f"{pair_id}.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {receptor}\n      msa: {a3m_rel}\n"
        f"  - protein:\n      id: B\n      sequence: {epitope}\n      msa: empty\n")
    rf3.append({"name": pair_id, "components": [
        {"seq": receptor, "chain_id": "A", "msa_path": a3m_rel},
        {"seq": epitope, "chain_id": "B"}]})

    fa_path = cf_fa_dir / f"{pair_id}.fasta"
    fa_path.write_text(f">{pair_id}\n{receptor}:{epitope}\n")
    cf_a3m_path = cf_msa_dir / f"{pair_id}.a3m"
    cf_a3m_path.write_text(complex_unpaired_a3m(receptor, epitope, a3m_text))

    manifest.append([pair_id, row["row_id"], locus, row["Ligand"], a3m_rel,
                     portable(fa_path), portable(cf_a3m_path), "built"])

(out / "rf3.json").write_text(json.dumps(rf3, indent=1))
csv.writer(open(out / "manifest.tsv", "w"), delimiter="\t").writerows(manifest)

skipped = [m for m in manifest[1:] if m[-1] != "built"]
print(f"built {len(rf3)}, skipped {len(skipped)} -> {portable(out)}")
print(f"  boltz/          {len(rf3)} yaml")
print(f"  colabfold/      {len(rf3)} fasta  (receptor:peptide)")
print(f"  colabfold_msa/  {len(rf3)} a3m    (unpaired custom MSA; preferred for AF2)")
print(f"  rf3.json")
for m in skipped:
    print("  ", m[2], m[-1])
