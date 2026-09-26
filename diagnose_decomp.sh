#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# DIAGNOSTIC ONLY - reads files, changes nothing.
#
# "DecompError: Mismatch in number of decomp terms!" means the
# complex, receptor and ligand decompositions did not produce
# the same set of residues. This gathers the facts needed to
# say exactly why, instead of guessing.
# ============================================================

set -Eeuo pipefail

SYS="${1:-mmgbsa_checks_focused/5sa23_HIE}"

echo "============================================================"
echo "DECOMPOSITION DIAGNOSTIC"
echo "System: $SYS"
echo "============================================================"

[[ -d "$SYS" ]] || { echo "ERROR: $SYS not found."; exit 1; }

echo
echo "--- 1. print_res actually used ---"
if [[ -s "$SYS/mmgbsa.in" ]]; then
    cat "$SYS/mmgbsa.in"
else
    echo "mmgbsa.in not found (it may have been overwritten by the retry)."
fi

echo
echo "--- 2. residue counts in each topology ---"
python3 - "$SYS" <<'PY'
import re
import sys
from pathlib import Path

sysdir = Path(sys.argv[1])


def residues(prmtop: Path):
    """Residue labels and their first-atom pointers from an Amber prmtop."""
    if not prmtop.is_file():
        return None, None
    text = prmtop.read_text(errors="replace")

    def flag(name):
        m = re.search(
            rf"%FLAG {name}\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG|\Z)",
            text, re.S)
        return m.group(1) if m else ""

    labels = flag("RESIDUE_LABEL").split()
    pointers = flag("RESIDUE_POINTER").split()
    return labels, pointers


for name in ("complex", "receptor", "ligand"):
    labels, _ = residues(sysdir / f"{name}.prmtop")
    if labels is None:
        print(f"  {name:<9} prmtop missing")
        continue
    print(f"  {name:<9} {len(labels):>5} residues   last: {' '.join(labels[-3:])}")

labels, _ = residues(sysdir / "complex.prmtop")
if labels:
    # Where do FMN and LIG sit in the complex numbering?
    for target in ("FMN", "LIG"):
        hits = [i + 1 for i, l in enumerate(labels) if l == target]
        print(f"  {target} residue number(s) in complex: "
              f"{hits if hits else 'NOT FOUND'}")

rec, _ = residues(sysdir / "receptor.prmtop")
lig, _ = residues(sysdir / "ligand.prmtop")
if labels and rec and lig:
    print(f"  complex = receptor + ligand ? "
          f"{len(labels)} vs {len(rec)} + {len(lig)} = {len(rec) + len(lig)}"
          f"  -> {'consistent' if len(labels) == len(rec) + len(lig) else 'MISMATCH'}")
PY

echo
echo "--- 3. residues requested, and where they fall ---"
python3 - "$SYS" <<'PY'
import re
import sys
from pathlib import Path

sysdir = Path(sys.argv[1])
inp = sysdir / "mmgbsa.in"
if not inp.is_file():
    raise SystemExit("  mmgbsa.in not available")

m = re.search(r'print_res\s*=\s*"([^"]*)"', inp.read_text())
if not m:
    print("  no print_res found (decomposition may have been dropped on retry)")
    raise SystemExit(0)

spec = m.group(1)
nums = set()
for part in spec.split(","):
    part = part.strip()
    if "-" in part:
        a, b = part.split("-")
        nums.update(range(int(a), int(b) + 1))
    elif part.isdigit():
        nums.add(int(part))

text = (sysdir / "complex.prmtop").read_text(errors="replace")
mm = re.search(r"%FLAG RESIDUE_LABEL\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG|\Z)",
               text, re.S)
labels = mm.group(1).split() if mm else []

rec_text = (sysdir / "receptor.prmtop").read_text(errors="replace")
mr = re.search(r"%FLAG RESIDUE_LABEL\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG|\Z)",
               rec_text, re.S)
n_rec = len(mr.group(1).split()) if mr else 0

print(f"  print_res = \"{spec}\"")
print(f"  residues requested: {len(nums)}")
print(f"  receptor spans residues 1-{n_rec}; ligand is residue {n_rec + 1}"
      f" of {len(labels)}")

in_rec = sorted(n for n in nums if n <= n_rec)
in_lig = sorted(n for n in nums if n > n_rec)
print(f"    requested inside receptor: {len(in_rec)}")
print(f"    requested inside ligand  : {len(in_lig)} {in_lig if in_lig else ''}")

if not in_lig:
    print()
    print("  >> The ligand residue was NOT requested.")
    print("     MMPBSA.py decomposes the complex over these residues, but the")
    print("     ligand decomposition has nothing to match, so the term counts")
    print("     differ and it raises DecompError.")
else:
    print()
    print("  >> The ligand WAS included, so the mismatch has another cause.")

missing = [n for n in nums if n > len(labels)]
if missing:
    print(f"  >> Requested residues beyond the end of the complex: {missing}")
PY

echo
echo "--- 4. what MMPBSA.py reported ---"
if [[ -s "$SYS/mmpbsa.log" ]]; then
    grep -Ei 'decomp|mismatch|error|warning|print_res|residue' "$SYS/mmpbsa.log" \
        | head -n 20
fi

echo
echo "--- 5. retained intermediate files ---"
ls "$SYS" | grep -i '_MMPBSA_' | head -n 15 || echo "  none retained"

echo
echo "--- 6. AmberTools version ---"
MMPBSA.py --version 2>&1 | head -n 3 || true

echo
echo "============================================================"
echo "Nothing was changed. This was read-only."
echo "============================================================"
