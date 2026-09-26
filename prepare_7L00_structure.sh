#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 1: Download PDB 7L00 (C. difficile FabK + FMN + XCJ)
#         and separate protein / FMN / XCJ / water / other.
#         No PyMOL: bash + Python only.
# ============================================================

set -Eeuo pipefail

PDB_ID="7L00"
TXT_DIR="txt_outputs"
mkdir -p "$TXT_DIR"

header_txt() {
    echo "============================================================"
    echo "siam chowdhury"
    echo "Computational and Medicinal Chemistry"
    echo "[Dr. Alam's Research Team] www.alamresearch.org"
    echo "Arkansas State University"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"
    echo
}

echo "============================================================"
echo "7L00 STRUCTURE DOWNLOAD AND SEPARATION"
echo "============================================================"

# ------------------------------------------------------------
# 1. Download (skipped if files already exist)
# ------------------------------------------------------------

download() {
    local url="$1"
    local out="$2"

    if [[ -s "$out" ]]; then
        echo "Found existing $out (download skipped)"
        return
    fi

    if command -v wget >/dev/null 2>&1; then
        wget -q -O "$out" "$url"
    else
        curl -fsSL -o "$out" "$url"
    fi
}

download "https://files.rcsb.org/download/${PDB_ID}.pdb" "${PDB_ID}.pdb"
download "https://files.rcsb.org/download/${PDB_ID}.cif" "${PDB_ID}.cif" \
    || echo "NOTE: ${PDB_ID}.cif not downloaded (optional, not needed)."

if [[ ! -s "${PDB_ID}.pdb" ]]; then
    echo "ERROR: ${PDB_ID}.pdb missing or empty."
    exit 1
fi
[[ -s "${PDB_ID}.cif" ]] || rm -f "${PDB_ID}.cif"

{
    header_txt
    echo "01 DOWNLOAD REPORT: ${PDB_ID}"
    echo "--------------------------"
    ls -lh "${PDB_ID}.pdb" "${PDB_ID}.cif" 2>/dev/null || ls -lh "${PDB_ID}.pdb"
    echo
    echo "TITLE / RESOLUTION / SOURCE"
    echo "---------------------------"
    grep -E '^(TITLE|COMPND   2|SOURCE   2|REMARK   2 RESOLUTION)' \
        "${PDB_ID}.pdb" || true
} > "$TXT_DIR/01_${PDB_ID}_download_report.txt"

cat "$TXT_DIR/01_${PDB_ID}_download_report.txt"

# ------------------------------------------------------------
# 2. Separate components and write reports (Python)
# ------------------------------------------------------------

python3 - "$PDB_ID" "$TXT_DIR" <<'PY'
import math
import sys
from collections import OrderedDict, defaultdict
from datetime import datetime
from pathlib import Path

pdb_id = sys.argv[1]
txt_dir = Path(sys.argv[2])

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
STAMP = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "MSE", "HID", "HIE", "HIP", "CYX",
}
WATER = {"HOH", "WAT", "DOD"}
COFACTOR = "FMN"
LIGAND = "XCJ"


def txt_header():
    return ["=" * 60, *IDENTITY, f"Date: {STAMP}", "=" * 60, ""]


def pdb_header(title):
    lines = [f"REMARK 999 {line}" for line in IDENTITY]
    lines.append(f"REMARK 999 {title}")
    lines.append(f"REMARK 999 Created: {STAMP}")
    return lines


def xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def res_key(line):
    return (line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip())


# ---------- read first model, resolve altlocs ----------
records = []
altloc_dropped = 0

for line in Path(f"{pdb_id}.pdb").read_text(errors="replace").splitlines():
    if line.startswith("ENDMDL"):
        break
    if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
        continue
    alt = line[16]
    if alt not in (" ", "A"):
        altloc_dropped += 1
        continue
    records.append(line[:16] + " " + line[17:])

if not records:
    raise SystemExit("ERROR: no coordinate records found.")

groups = {"protein": [], "fmn": [], "xcj": [], "water": [], "other": []}

for line in records:
    name = line[17:20].strip()
    if name in AMINO:
        groups["protein"].append(line)
    elif name == COFACTOR:
        groups["fmn"].append(line)
    elif name == LIGAND:
        groups["xcj"].append(line)
    elif name in WATER:
        groups["water"].append(line)
    else:
        groups["other"].append(line)


def write_pdb(path, lines, title, ter_by_chain=False):
    out = pdb_header(title)
    previous_chain = None
    for line in lines:
        if ter_by_chain and previous_chain is not None and line[21] != previous_chain:
            out.append("TER")
        out.append(line.rstrip())
        previous_chain = line[21]
    if lines:
        out.append("TER")
    out.append("END")
    Path(path).write_text("\n".join(out) + "\n")


write_pdb(f"{pdb_id}_protein.pdb", groups["protein"], f"{pdb_id} protein only", True)
write_pdb(f"{pdb_id}_FMN.pdb", groups["fmn"], f"{pdb_id} FMN cofactor", True)
write_pdb(f"{pdb_id}_XCJ_all.pdb", groups["xcj"], f"{pdb_id} native XCJ (all copies)", True)
write_pdb(f"{pdb_id}_water.pdb", groups["water"], f"{pdb_id} crystal waters")
write_pdb(f"{pdb_id}_other.pdb", groups["other"], f"{pdb_id} other components")
write_pdb(
    f"{pdb_id}_protein_FMN.pdb",
    groups["protein"] + groups["fmn"],
    f"{pdb_id} receptor draft: protein + FMN (XCJ, water removed)",
    True,
)

# ---------- 02 component inventory ----------
het_counts = defaultdict(int)
het_residues = defaultdict(set)
for line in records:
    if line.startswith("HETATM"):
        het_counts[line[17:20].strip()] += 1
        het_residues[line[17:20].strip()].add(res_key(line))

report = txt_header() + [
    f"02 COMPONENT INVENTORY: {pdb_id} (model 1, altloc A/blank kept)",
    "-" * 60,
    f"Total coordinate records kept : {len(records)}",
    f"Alternate-location atoms dropped: {altloc_dropped}",
    "",
    f"{'GROUP':<10}{'ATOMS':>8}   OUTPUT FILE",
    f"{'protein':<10}{len(groups['protein']):>8}   {pdb_id}_protein.pdb",
    f"{'FMN':<10}{len(groups['fmn']):>8}   {pdb_id}_FMN.pdb",
    f"{'XCJ':<10}{len(groups['xcj']):>8}   {pdb_id}_XCJ_all.pdb",
    f"{'water':<10}{len(groups['water']):>8}   {pdb_id}_water.pdb",
    f"{'other':<10}{len(groups['other']):>8}   {pdb_id}_other.pdb",
    "",
    f"Receptor draft (protein + FMN): {pdb_id}_protein_FMN.pdb "
    f"({len(groups['protein']) + len(groups['fmn'])} atoms)",
    "",
    "HETATM COMPONENTS",
    f"{'RESNAME':<10}{'ATOMS':>8}{'COPIES':>8}",
]
for name in sorted(het_counts):
    report.append(f"{name:<10}{het_counts[name]:>8}{len(het_residues[name]):>8}")

if not groups["xcj"]:
    report += ["", "WARNING: no XCJ found. Check the HETATM list above."]
if not groups["fmn"]:
    report += ["", "WARNING: no FMN found. Check the HETATM list above."]

(txt_dir / f"02_{pdb_id}_component_inventory.txt").write_text("\n".join(report) + "\n")

# ---------- 03 chain report ----------
chain_atoms = defaultdict(int)
chain_res = defaultdict(OrderedDict)
for line in groups["protein"]:
    chain_atoms[line[21]] += 1
    chain_res[line[21]][res_key(line)] = True

report = txt_header() + [
    f"03 PROTEIN CHAIN REPORT: {pdb_id}",
    "-" * 60,
    f"{'CHAIN':<8}{'RESIDUES':>10}{'ATOMS':>10}{'FIRST':>8}{'LAST':>8}",
]
for ch in chain_atoms:
    keys = list(chain_res[ch])
    report.append(
        f"{ch:<8}{len(keys):>10}{chain_atoms[ch]:>10}"
        f"{keys[0][1]:>8}{keys[-1][1]:>8}"
    )
(txt_dir / f"03_{pdb_id}_chain_report.txt").write_text("\n".join(report) + "\n")


# ---------- 04/05 per-copy geometry for XCJ and FMN ----------
def copies(lines):
    out = OrderedDict()
    for line in lines:
        out.setdefault(res_key(line), []).append(line)
    return out


def geom(lines):
    pts = [xyz(line) for line in lines]
    xs, ys, zs = zip(*pts)
    sx, sy, sz = max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)
    return {
        "n": len(pts),
        "c": (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)),
        "min": (min(xs), min(ys), min(zs)),
        "max": (max(xs), max(ys), max(zs)),
        "size": (sx, sy, sz),
        "diag": math.sqrt(sx * sx + sy * sy + sz * sz),
    }


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


xcj_copies = copies(groups["xcj"])
fmn_copies = copies(groups["fmn"])
protein_xyz = [(line[21], xyz(line)) for line in groups["protein"]]


def copy_report(title, copy_map, file_tag):
    lines = txt_header() + [title, "-" * 60, ""]
    for key, atoms in copy_map.items():
        g = geom(atoms)
        ch, num, ins, name = key

        # nearest protein chain (minimum atom-atom distance)
        near = defaultdict(lambda: float("inf"))
        for pch, p in protein_xyz:
            for a in atoms:
                d = dist(p, xyz(a))
                if d < near[pch]:
                    near[pch] = d

        lines += [
            f"RESNAME {name}   CHAIN {ch}   RESIDUE {num}{ins}",
            f"ATOM_COUNT {g['n']}",
            f"CENTER_X {g['c'][0]:.3f}",
            f"CENTER_Y {g['c'][1]:.3f}",
            f"CENTER_Z {g['c'][2]:.3f}",
            f"X_MIN {g['min'][0]:.3f}   X_MAX {g['max'][0]:.3f}",
            f"Y_MIN {g['min'][1]:.3f}   Y_MAX {g['max'][1]:.3f}",
            f"Z_MIN {g['min'][2]:.3f}   Z_MAX {g['max'][2]:.3f}",
            f"SIZE_X {g['size'][0]:.3f} A",
            f"SIZE_Y {g['size'][1]:.3f} A",
            f"SIZE_Z {g['size'][2]:.3f} A",
            f"BOUNDING_BOX_DIAGONAL {g['diag']:.3f} A",
            "NEAREST PROTEIN CHAINS (min atom distance):",
        ]
        for pch, d in sorted(near.items(), key=lambda item: item[1])[:3]:
            lines.append(f"   chain {pch}: {d:.3f} A")
        lines.append("")

        tag = f"{name}_chain{ch}_res{num}{ins}"
        write_pdb(f"{pdb_id}_{tag}.pdb", atoms, f"{pdb_id} {name} chain {ch} residue {num}{ins}")
        lines.append(f"FILE {pdb_id}_{tag}.pdb")
        lines.append("")

    (txt_dir / file_tag).write_text("\n".join(lines) + "\n")


copy_report(f"04 XCJ NATIVE-LIGAND COPIES: {pdb_id}", xcj_copies, f"04_{pdb_id}_XCJ_copies.txt")
copy_report(f"05 FMN COFACTOR COPIES: {pdb_id}", fmn_copies, f"05_{pdb_id}_FMN_copies.txt")

# ---------- 06 XCJ-FMN pairing (which FMN belongs to which pocket) ----------
lines = txt_header() + [
    f"06 XCJ-FMN PAIRING: {pdb_id}",
    "-" * 60,
    "Minimum atom-atom distance between each XCJ copy and each FMN copy.",
    "The closest FMN defines the active site used for that XCJ copy.",
    "",
    f"{'XCJ':<18}{'FMN':<18}{'MIN_DIST_A':>12}{'CENTER_DIST_A':>15}",
]
for xk, xa in xcj_copies.items():
    for fk, fa in fmn_copies.items():
        dmin = min(dist(xyz(a), xyz(b)) for a in xa for b in fa)
        dcen = dist(geom(xa)["c"], geom(fa)["c"])
        lines.append(
            f"{xk[0] + ':' + xk[1]:<18}{fk[0] + ':' + fk[1]:<18}"
            f"{dmin:>12.3f}{dcen:>15.3f}"
        )
(txt_dir / f"06_{pdb_id}_XCJ_FMN_pairing.txt").write_text("\n".join(lines) + "\n")

print(f"Protein atoms : {len(groups['protein'])}")
print(f"FMN atoms     : {len(groups['fmn'])} in {len(fmn_copies)} copies")
print(f"XCJ atoms     : {len(groups['xcj'])} in {len(xcj_copies)} copies")
print(f"Water atoms   : {len(groups['water'])}")
print(f"Other atoms   : {len(groups['other'])}")
PY

echo
echo "============================================================"
echo "STEP 1 COMPLETED. Read:"
echo "============================================================"
echo "cat $TXT_DIR/02_${PDB_ID}_component_inventory.txt"
echo "cat $TXT_DIR/03_${PDB_ID}_chain_report.txt"
echo "cat $TXT_DIR/04_${PDB_ID}_XCJ_copies.txt"
echo "cat $TXT_DIR/05_${PDB_ID}_FMN_copies.txt"
echo "cat $TXT_DIR/06_${PDB_ID}_XCJ_FMN_pairing.txt"
