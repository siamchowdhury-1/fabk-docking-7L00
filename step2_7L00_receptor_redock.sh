#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 2: 7L00 FabK receptor preparation (chains A+B + FMN),
#         active-site characterization, XCJ redocking and
#         RMSD validation of the docking protocol.
# ============================================================

set -Eeuo pipefail

# ------------------------- settings -------------------------
PDB_ID="7L00"
SITE_CHAIN="A"          # pocket chain (XCJ A:402, FMN A:401)
PARTNER_CHAIN="B"       # dimer partner lining the pocket
XCJ_RESSEQ="402"
FMN_RESSEQ="401"
PH="7.4"
BOX_PADDING="10.0"      # added to largest XCJ span (Angstrom)
EXHAUSTIVENESS="32"
NUM_MODES="20"
ENERGY_RANGE="5"
SEED="12345"
SITE_CUTOFF="5.0"       # residue cutoff for site characterization
TXT_DIR="txt_outputs"
# ------------------------------------------------------------

mkdir -p "$TXT_DIR"
TAG="${PDB_ID}_${SITE_CHAIN}${PARTNER_CHAIN}"

IDENTITY=(
    "siam chowdhury"
    "Computational and Medicinal Chemistry"
    "[Dr. Alam's Research Team] www.alamresearch.org"
    "Arkansas State University"
)

header_txt() {
    echo "============================================================"
    printf '%s\n' "${IDENTITY[@]}"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"
    echo
}

stamp_pdb() {   # REMARK 999 lines at top (PDB / PDBQT)
    local f="$1" title="$2" tmp
    tmp="$(mktemp)"
    {
        for l in "${IDENTITY[@]}"; do echo "REMARK 999 $l"; done
        echo "REMARK 999 $title"
        echo "REMARK 999 Created: $(date '+%Y-%m-%d %H:%M:%S')"
        grep -v '^REMARK 999 ' "$f" || true
    } > "$tmp"
    mv "$tmp" "$f"
}

stamp_mol2() {  # '#' comment lines at top
    local f="$1" title="$2" tmp
    tmp="$(mktemp)"
    {
        for l in "${IDENTITY[@]}"; do echo "# $l"; done
        echo "# $title"
        grep -v '^# ' "$f" || true
    } > "$tmp"
    mv "$tmp" "$f"
}

stamp_sdf() {   # PREPARED_BY data field in every record
    local f="$1" tmp
    tmp="$(mktemp)"
    awk -v a="${IDENTITY[0]}" -v b="${IDENTITY[1]}" \
        -v c="${IDENTITY[2]}" -v d="${IDENTITY[3]}" '
        /^> *<PREPARED_BY>/ {skip=1; next}
        skip && /^$/ {skip=0; next}
        skip {next}
        /^\$\$\$\$/ {
            print "> <PREPARED_BY>"
            print a; print b; print c; print d
            print ""
        }
        {print}
    ' "$f" > "$tmp"
    mv "$tmp" "$f"
}

echo "============================================================"
echo "STEP 2: ${PDB_ID} RECEPTOR + XCJ REDOCKING VALIDATION"
echo "============================================================"

for tool in python3 obabel vina mk_prepare_ligand.py; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not found. Run: conda activate docking"
        exit 1
    }
done

for f in "${PDB_ID}_protein.pdb" "${PDB_ID}_FMN.pdb" \
         "${PDB_ID}_XCJ_chain${SITE_CHAIN}_res${XCJ_RESSEQ}.pdb"; do
    [[ -s "$f" ]] || { echo "ERROR: $f missing. Run Step 1 first."; exit 1; }
done

# ============================================================
# 2A. Build receptor PDB files (chains A+B, FMN A+B)
# ============================================================

python3 - "$PDB_ID" "$SITE_CHAIN" "$PARTNER_CHAIN" "$TAG" <<'PY'
import sys
from pathlib import Path

pdb_id, site, partner, tag = sys.argv[1:5]
chains = {site, partner}


def pick(path):
    return [
        l for l in Path(path).read_text().splitlines()
        if l.startswith(("ATOM", "HETATM")) and l[21] in chains
    ]


def write(path, lines):
    out, prev = [], None
    for l in lines:
        if prev is not None and l[21] != prev:
            out.append("TER")
        out.append(l)
        prev = l[21]
    out += ["TER", "END"]
    Path(path).write_text("\n".join(out) + "\n")


protein = pick(f"{pdb_id}_protein.pdb")
fmn = pick(f"{pdb_id}_FMN.pdb")

# Residues with missing side-chain atoms have no Meeko template and
# would be deleted. Trim them to a CB stub (ALA) or backbone (GLY).
EXPECTED = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9,
    "GLU": 9, "GLY": 4, "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9,
    "MET": 8, "PHE": 11, "PRO": 7, "SER": 6, "THR": 7, "TRP": 14,
    "TYR": 12, "VAL": 7,
}
groups = {}
for l in protein:
    groups.setdefault((l[21], l[22:27]), []).append(l)

xcj = [
    (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    for l in Path(f"{pdb_id}_XCJ_chain{site}_res402.pdb").read_text().splitlines()
    if l.startswith(("ATOM", "HETATM"))
]

fixed_lines, notes = [], []
for (ch, num), lines in groups.items():
    name = lines[0][17:20]
    atoms = {l[12:16].strip() for l in lines if l[12:16].strip() != "OXT"}
    need = EXPECTED.get(name)
    if need is not None and len(atoms) < need and name != "PRO":
        keep = {"N", "CA", "C", "O", "OXT", "CB"}
        new = "ALA" if "CB" in atoms else "GLY"
        lines = [l[:17] + new + l[20:] for l in lines if l[12:16].strip() in keep]
        dmin = min(
            ((float(l[30:38]) - x) ** 2 + (float(l[38:46]) - y) ** 2
             + (float(l[46:54]) - z) ** 2) ** 0.5
            for l in lines for (x, y, z) in xcj
        )
        notes.append(
            f"{ch}:{num.strip():<6}{name} -> {new}  "
            f"(had {len(atoms)}/{need} heavy atoms, {dmin:.1f} A from XCJ {site})"
        )
    fixed_lines += lines
protein = fixed_lines
Path(f"{tag}_truncated_residues.txt").write_text(
    "\n".join(notes) + ("\n" if notes else "None\n")
)

write(f"{tag}_protein.pdb", protein)
write(f"{tag}_FMN.pdb", fmn)
write(f"{tag}_protein_FMN.pdb", protein + fmn)

print(f"Protein atoms (chains {''.join(sorted(chains))}): {len(protein)}")
print(f"FMN atoms: {len(fmn)}")
PY

{
    header_txt
    echo "07 RECEPTOR BUILD: ${TAG}"
    echo "------------------------------------------------------------"
    echo "Pocket chain         : ${SITE_CHAIN} (XCJ ${SITE_CHAIN}:${XCJ_RESSEQ}, FMN ${SITE_CHAIN}:${FMN_RESSEQ})"
    echo "Partner chain        : ${PARTNER_CHAIN} (XCJ-${PARTNER_CHAIN} min distance ~7 A)"
    echo "Kept                 : protein chains ${SITE_CHAIN}+${PARTNER_CHAIN}, FMN of both chains"
    echo "Removed              : all XCJ, waters (HOH), GOL, NA, chains C and D"
    echo
    echo "Incomplete side chains trimmed (crystal atoms missing):"
    sed 's/^/    /' "${TAG}_truncated_residues.txt"
    echo
    for f in "${TAG}_protein.pdb" "${TAG}_FMN.pdb" "${TAG}_protein_FMN.pdb"; do
        printf '%-32s %6s atoms\n' "$f" "$(grep -Ec '^(ATOM|HETATM)' "$f")"
    done
} > "$TXT_DIR/07_${TAG}_receptor_build.txt"
cat "$TXT_DIR/07_${TAG}_receptor_build.txt"

# ============================================================
# 2B. Crystal reference ligands with correct bond orders
#     (RCSB ModelServer SDF; fallback: Open Babel perception)
# ============================================================

sdf_heavy() {   # heavy-atom count of the first SDF record
    awk 'NR==4 {n=substr($0,1,3)+0} NR>4 && NR<=4+n && $4!="H" {h++}
         NR>4+n {exit} END {print h+0}' "$1"
}

fetch_ligand_sdf() {   # chain resseq expected_heavy out_sdf fallback_pdb
    local ch="$1" num="$2" expect="$3" out="$4" fallback="$5"
    local url="https://models.rcsb.org/v1/${PDB_ID,,}/ligand?auth_asym_id=${ch}&auth_seq_id=${num}&encoding=sdf"
    local source="RCSB ModelServer (CCD bond orders)"

    rm -f "$out"
    if command -v wget >/dev/null 2>&1; then
        wget -q -O "$out" "$url" || true
    else
        curl -fsSL -o "$out" "$url" || true
    fi

    local heavy=0
    [[ -s "$out" ]] && heavy=$(sdf_heavy "$out")

    if [[ "${heavy:-0}" -ne "$expect" ]]; then
        source="Open Babel bond perception from crystal PDB (fallback)"
        obabel "$fallback" -O "$out" 2>/dev/null
        heavy=$(sdf_heavy "$out")
    fi

    echo "$source|$heavy"
}

XCJ_PDB="${PDB_ID}_XCJ_chain${SITE_CHAIN}_res${XCJ_RESSEQ}.pdb"
XCJ_REF="XCJ_${SITE_CHAIN}${XCJ_RESSEQ}_crystal.sdf"
XCJ_INFO=$(fetch_ligand_sdf "$SITE_CHAIN" "$XCJ_RESSEQ" 27 "$XCJ_REF" "$XCJ_PDB")

{
    header_txt
    echo "12 XCJ CRYSTAL REFERENCE"
    echo "------------------------------------------------------------"
    echo "Copy         : XCJ ${SITE_CHAIN}:${XCJ_RESSEQ}"
    echo "Source       : ${XCJ_INFO%%|*}"
    echo "Heavy atoms  : ${XCJ_INFO##*|} (expected 27)"
    echo "File         : ${XCJ_REF}"
    echo
    echo "SMILES (crystal reference):"
    obabel "$XCJ_REF" -osmi 2>/dev/null
} > "$TXT_DIR/12_XCJ_crystal_reference.txt"
cat "$TXT_DIR/12_XCJ_crystal_reference.txt"

# ============================================================
# 2C. Protein PDBQT (Meeko; fallback Open Babel)
# ============================================================

PROT_BASE="${TAG}_protein_prep"
PROT_REPORT="$TXT_DIR/09_${TAG}_protein_prepare_report.txt"
PROT_METHOD="Meeko mk_prepare_receptor.py"

rm -f "${PROT_BASE}.pdbqt"
{ header_txt; echo "09 PROTEIN PDBQT PREPARATION"; echo; } > "$PROT_REPORT"

if command -v mk_prepare_receptor.py >/dev/null 2>&1; then
    HELP_OUTPUT="$(mk_prepare_receptor.py -h 2>&1 || true)"

    if grep -q -- "--read_pdb" <<< "$HELP_OUTPUT"; then
        in_args=(--read_pdb "${TAG}_protein.pdb")
    else
        in_args=(-i "${TAG}_protein.pdb")
    fi

    bad_args=()
    if grep -q -- "--allow_bad_res" <<< "$HELP_OUTPUT"; then
        bad_args=(--allow_bad_res)
    elif grep -q -- "--delete_bad_res" <<< "$HELP_OUTPUT"; then
        bad_args=(--delete_bad_res)
    fi

    echo "COMMAND: mk_prepare_receptor.py ${in_args[*]} -o ${PROT_BASE} -p ${bad_args[*]} --default_altloc A" \
        >> "$PROT_REPORT"
    echo >> "$PROT_REPORT"

    mk_prepare_receptor.py "${in_args[@]}" -o "$PROT_BASE" -p \
        "${bad_args[@]}" --default_altloc A \
        >> "$PROT_REPORT" 2>&1 || true
fi

if [[ ! -s "${PROT_BASE}.pdbqt" ]]; then
    PROT_METHOD="Open Babel (-xr, pH ${PH}, Gasteiger) fallback"
    echo "Meeko failed or unavailable -> Open Babel fallback" >> "$PROT_REPORT"
    obabel "${TAG}_protein.pdb" -O "${PROT_BASE}.pdbqt" \
        -xr -p "$PH" --partialcharge gasteiger >> "$PROT_REPORT" 2>&1
fi

[[ -s "${PROT_BASE}.pdbqt" ]] || { echo "ERROR: protein PDBQT not created."; exit 1; }
echo "METHOD USED: $PROT_METHOD" >> "$PROT_REPORT"

# ============================================================
# 2D. FMN PDBQT (bond orders from RCSB, H at pH 7.4)
# ============================================================

FMN_REPORT="$TXT_DIR/10_${TAG}_FMN_prepare_report.txt"
{ header_txt; echo "10 FMN COFACTOR PREPARATION (kept in receptor)"; echo; } > "$FMN_REPORT"

rm -f "${TAG}_FMN_combined.pdbqt"

for ch in "$SITE_CHAIN" "$PARTNER_CHAIN"; do
    fmn_pdb="FMN_${ch}${FMN_RESSEQ}_crystal.pdb"
    fmn_sdf="FMN_${ch}${FMN_RESSEQ}_crystal.sdf"
    fmn_qt="FMN_${ch}${FMN_RESSEQ}.pdbqt"

    grep -E '^(ATOM|HETATM)' "${PDB_ID}_FMN.pdb" | awk -v c="$ch" 'substr($0,22,1)==c' > "$fmn_pdb"
    echo "END" >> "$fmn_pdb"

    info=$(fetch_ligand_sdf "$ch" "$FMN_RESSEQ" 31 "$fmn_sdf" "$fmn_pdb")

    obabel "$fmn_sdf" -O "$fmn_qt" -xr -p "$PH" \
        --partialcharge gasteiger >> "$FMN_REPORT" 2>&1

    {
        echo "FMN ${ch}:${FMN_RESSEQ}"
        echo "  source      : ${info%%|*}"
        echo "  heavy atoms : ${info##*|} (expected 31)"
        echo "  PDBQT atoms : $(grep -Ec '^(ATOM|HETATM)' "$fmn_qt")"
        echo "  SMILES      : $(obabel "$fmn_sdf" -osmi 2>/dev/null | awk '{print $1}')"
        echo
    } >> "$FMN_REPORT"

    # rename residue to FMN / chain / resseq, keep coordinates + types
    awk -v c="$ch" -v r="$FMN_RESSEQ" '
        /^(ATOM|HETATM)/ {
            printf "HETATM%sFMN %s%4d%s\n",
                substr($0,7,11), c, r, substr($0,27)
        }' "$fmn_qt" >> "${TAG}_FMN_combined.pdbqt"
done

# ============================================================
# 2E. Final receptor PDBQT = protein + FMN
# ============================================================

RECEPTOR="${TAG}_receptor.pdbqt"
{
    grep -E '^(ATOM|HETATM)' "${PROT_BASE}.pdbqt"
    cat "${TAG}_FMN_combined.pdbqt"
} > "$RECEPTOR"
stamp_pdb "$RECEPTOR" "${PDB_ID} receptor: chains ${SITE_CHAIN}+${PARTNER_CHAIN} + FMN (protein: ${PROT_METHOD})"

for f in "${TAG}_protein.pdb" "${TAG}_FMN.pdb" "${TAG}_protein_FMN.pdb" \
         "${PROT_BASE}.pdbqt"; do
    stamp_pdb "$f" "${PDB_ID} receptor component"
done
for f in FMN_*_crystal.pdb FMN_*.pdbqt; do stamp_pdb "$f" "${PDB_ID} FMN cofactor"; done
for f in FMN_*_crystal.sdf; do stamp_sdf "$f"; done

{
    header_txt
    echo "11 FINAL RECEPTOR PDBQT CHECK: ${RECEPTOR}"
    echo "------------------------------------------------------------"
    ls -lh "$RECEPTOR"
    echo
    echo "Protein method       : ${PROT_METHOD}"
    echo "Input heavy atoms    : $(grep -Ec '^(ATOM|HETATM)' "${TAG}_protein_FMN.pdb")"
    echo "PDBQT atoms total    : $(grep -Ec '^(ATOM|HETATM)' "$RECEPTOR")"
    echo "FMN atoms in PDBQT   : $(grep -E '^(ATOM|HETATM)' "$RECEPTOR" | awk 'substr($0,18,3)=="FMN"' | wc -l)"
    echo
    echo "ATOM TYPES"
    awk '/^(ATOM|HETATM)/ {n[$NF]++} END {for (t in n) printf "%8d %s\n", n[t], t}' \
        "$RECEPTOR" | sort -k2
    echo
    echo "WARNINGS / DELETED RESIDUES (protein preparation)"
    grep -Ein 'error|warning|no template|unknown|cannot|deleted|skip' "$PROT_REPORT" \
        || echo "None detected."
} > "$TXT_DIR/11_${TAG}_receptor_PDBQT_check.txt"
cat "$TXT_DIR/11_${TAG}_receptor_PDBQT_check.txt"

# ============================================================
# 2F. Active-site characterization + grid box (from crystal XCJ)
# ============================================================

python3 - "$XCJ_PDB" "${TAG}_protein_FMN.pdb" "$SITE_CUTOFF" "$BOX_PADDING" \
          "$TXT_DIR" "$TAG" "$SITE_CHAIN" "$XCJ_RESSEQ" <<'PY'
import math
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

xcj_path, rec_path, cutoff, pad, txt_dir, tag, site, resseq = sys.argv[1:9]
cutoff, pad = float(cutoff), float(pad)
txt_dir = Path(txt_dir)

HEAD = [
    "=" * 60,
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
    f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}",
    "=" * 60,
    "",
]


def atoms(path):
    out = []
    for l in Path(path).read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            out.append({
                "name": l[12:16].strip(),
                "res": l[17:20].strip(),
                "chain": l[21],
                "num": l[22:26].strip(),
                "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
            })
    return out


lig = atoms(xcj_path)
rec = atoms(rec_path)
xs, ys, zs = zip(*[a["xyz"] for a in lig])
center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
span = (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
size = math.ceil(max(span) + pad)

Path(f"{tag}_grid_box.txt").write_text(
    f"CENTER_X {center[0]:.3f}\nCENTER_Y {center[1]:.3f}\nCENTER_Z {center[2]:.3f}\n"
    f"SIZE_X {size:.1f}\nSIZE_Y {size:.1f}\nSIZE_Z {size:.1f}\n"
)

report = HEAD + [
    f"14 GRID BOX DEFINITION: {tag}",
    "-" * 60,
    f"Reference ligand : crystal XCJ {site}:{resseq}",
    f"Center           : {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}",
    f"XCJ span         : {span[0]:.3f} x {span[1]:.3f} x {span[2]:.3f} A",
    f"Padding          : {pad:.1f} A added to the largest span",
    f"Box size         : {size} x {size} x {size} A",
    f"Box file         : {tag}_grid_box.txt",
]
(txt_dir / f"14_{tag}_grid_box.txt").write_text("\n".join(report) + "\n")

# residues within cutoff of any XCJ atom
res = OrderedDict()
for a in rec:
    d = min(math.dist(a["xyz"], b["xyz"]) for b in lig)
    key = (a["chain"], a["num"], a["res"])
    if key not in res or d < res[key][0]:
        res[key] = (d, a["name"])

site_res = sorted(
    [(k, v) for k, v in res.items() if v[0] <= cutoff],
    key=lambda item: item[1][0],
)

report = HEAD + [
    f"08 ACTIVE-SITE CHARACTERIZATION: {tag}",
    "-" * 60,
    f"Residues with any atom within {cutoff:.1f} A of crystal XCJ {site}:{resseq}",
    "(pocket used for all docking in this study)",
    "",
    f"{'CHAIN':<7}{'RESIDUE':<10}{'NUMBER':<8}{'MIN_DIST_A':>11}  CLOSEST_ATOM",
]
for (ch, num, name), (d, atom) in site_res:
    report.append(f"{ch:<7}{name:<10}{num:<8}{d:>11.3f}  {atom}")

by_chain = {}
for (ch, _, _), _ in site_res:
    by_chain[ch] = by_chain.get(ch, 0) + 1
report += ["", "RESIDUES PER CHAIN"] + [f"  chain {c}: {n}" for c, n in sorted(by_chain.items())]
if len(by_chain) > 1:
    report += [
        "",
        "Residues from more than one chain line the pocket (dimer interface).",
    ]
else:
    only = next(iter(by_chain)) if by_chain else "-"
    report += [
        "",
        f"All residues within {cutoff:.1f} A belong to chain {only}.",
        "The partner chain lies just outside this cutoff but inside the",
        "docking box, so it is kept to define the box boundary correctly.",
    ]
(txt_dir / f"08_{tag}_active_site_residues.txt").write_text("\n".join(report) + "\n")

print(f"Grid center {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}, size {size} A")
print(f"Active-site residues within {cutoff} A: {len(site_res)}")
PY

cat "$TXT_DIR/14_${TAG}_grid_box.txt"
cat "$TXT_DIR/08_${TAG}_active_site_residues.txt"

read -r CX CY CZ SX SY SZ < <(awk '{printf "%s ", $2} END {print ""}' "${TAG}_grid_box.txt")

# ============================================================
# 2G. XCJ ligand preparation (new 3D conformer: unbiased redock)
# ============================================================

XCJ_REPORT="$TXT_DIR/13_XCJ_ligand_prepare_report.txt"
{ header_txt; echo "13 XCJ LIGAND PREPARATION FOR REDOCKING"; echo; } > "$XCJ_REPORT"
echo "obabel ${XCJ_REF} -O XCJ_H.sdf -p ${PH} --gen3d --minimize --ff MMFF94 --steps 500" >> "$XCJ_REPORT"

obabel "$XCJ_REF" -O XCJ_H.sdf -p "$PH" --gen3d \
    --minimize --ff MMFF94 --steps 500 >> "$XCJ_REPORT" 2>&1

mk_prepare_ligand.py -i XCJ_H.sdf -o XCJ.pdbqt >> "$XCJ_REPORT" 2>&1
[[ -s XCJ.pdbqt ]] || { echo "ERROR: XCJ.pdbqt not created."; exit 1; }

{
    echo
    echo "SMILES after protonation (pH ${PH}):"
    obabel XCJ_H.sdf -osmi 2>/dev/null
    echo "PDBQT atoms : $(grep -Ec '^(ATOM|HETATM)' XCJ.pdbqt)"
    grep '^TORSDOF' XCJ.pdbqt || true
} >> "$XCJ_REPORT"
cat "$XCJ_REPORT"

stamp_sdf "$XCJ_REF"
stamp_sdf XCJ_H.sdf
stamp_pdb XCJ.pdbqt "XCJ ligand prepared for redocking"

# ============================================================
# 2H. Redocking with Vina
# ============================================================

CONFIG="${TAG}_XCJ_redock_config.txt"
OUT="XCJ_redock_out.pdbqt"

{
    for l in "${IDENTITY[@]}"; do echo "# $l"; done
    echo "# ${PDB_ID} XCJ redocking configuration"
    echo
    echo "receptor = ${RECEPTOR}"
    echo "ligand = XCJ.pdbqt"
    echo
    echo "center_x = ${CX}"
    echo "center_y = ${CY}"
    echo "center_z = ${CZ}"
    echo
    echo "size_x = ${SX}"
    echo "size_y = ${SY}"
    echo "size_z = ${SZ}"
    echo
    echo "exhaustiveness = ${EXHAUSTIVENESS}"
    echo "num_modes = ${NUM_MODES}"
    echo "energy_range = ${ENERGY_RANGE}"
} > "$CONFIG"

echo
echo "------------------------------------------------------------"
echo "REDOCKING XCJ INTO ${RECEPTOR}"
echo "------------------------------------------------------------"

{ header_txt; echo "15 VINA REDOCKING RUN (seed ${SEED})"; echo; } > "$TXT_DIR/15_XCJ_redock_vina_run.txt"
vina --config "$CONFIG" --out "$OUT" --seed "$SEED" 2>&1 \
    | tee -a "$TXT_DIR/15_XCJ_redock_vina_run.txt"

[[ -s "$OUT" ]] || { echo "ERROR: $OUT not created."; exit 1; }

{
    header_txt
    echo "16 XCJ REDOCKING SCORES"
    echo
    awk '/^mode[[:space:]]*\|/ {show=1} show {print}' "$TXT_DIR/15_XCJ_redock_vina_run.txt"
} > "$TXT_DIR/16_XCJ_redock_scores.txt"

# ============================================================
# 2I. Export poses (SDF, MOL2, PDB)
# ============================================================

rm -f XCJ_redock_poses.sdf
EXPORT_METHOD="Meeko mk_export.py"
if command -v mk_export.py >/dev/null 2>&1; then
    mk_export.py "$OUT" -s XCJ_redock_poses.sdf > /dev/null 2>&1 || true
fi
if [[ ! -s XCJ_redock_poses.sdf ]]; then
    EXPORT_METHOD="Open Babel (fallback)"
    obabel "$OUT" -O XCJ_redock_poses.sdf 2>/dev/null
fi

obabel XCJ_redock_poses.sdf -l 1 -O XCJ_redock_best.mol2 2>/dev/null
obabel XCJ_redock_poses.sdf -l 1 -O XCJ_redock_best.pdb 2>/dev/null
obabel "$XCJ_REF" -O XCJ_crystal.mol2 2>/dev/null

# ============================================================
# 2J. RMSD validation (symmetry-aware, heavy atoms, no fitting)
# ============================================================

python3 - "$XCJ_REF" XCJ_redock_poses.sdf "$TXT_DIR" "$EXPORT_METHOD" <<'PY'
import math
import sys
from datetime import datetime
from pathlib import Path

ref_path, poses_path, txt_dir, export_method = sys.argv[1:5]
txt_dir = Path(txt_dir)

HEAD = [
    "=" * 60,
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
    f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}",
    "=" * 60,
    "",
]

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")


def heavy_graph(mol):
    """Connectivity-only heavy-atom copy (ignores bond order, charge, H)."""
    conf = mol.GetConformer()
    keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    index = {old: new for new, old in enumerate(keep)}
    g = Chem.RWMol()
    for old in keep:
        g.AddAtom(Chem.Atom(mol.GetAtomWithIdx(old).GetAtomicNum()))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in index and j in index:
            g.AddBond(index[i], index[j], Chem.BondType.SINGLE)
    g = g.GetMol()
    g.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(g)
    xyz = [tuple(conf.GetAtomPosition(old)) for old in keep]
    return g, xyz


def load(path):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    return [m for m in supplier if m is not None]


ref = load(ref_path)[0]
ref_g, ref_xyz = heavy_graph(ref)
poses = load(poses_path)

rows = []
for n, pose in enumerate(poses, start=1):
    g, xyz = heavy_graph(pose)
    matches = g.GetSubstructMatches(ref_g, uniquify=False, maxMatches=100000)
    if not matches:
        rows.append((n, None, None))
        continue
    best = min(
        math.sqrt(sum(math.dist(ref_xyz[i], xyz[m[i]]) ** 2
                      for i in range(len(ref_xyz))) / len(ref_xyz))
        for m in matches
    )
    rows.append((n, best, len(matches)))

valid = [r for r in rows if r[1] is not None]
if not valid:
    raise SystemExit("ERROR: could not map docked poses onto crystal XCJ.")

top = rows[0][1]
best_row = min(valid, key=lambda r: r[1])

if top is None:
    verdict = "FAILED: pose 1 could not be matched"
elif top < 2.0:
    verdict = "EXCELLENT (< 2.0 A) - docking protocol VALIDATED"
elif top <= 3.0:
    verdict = "ACCEPTABLE (2-3 A) - protocol usable"
else:
    verdict = "RE-OPTIMIZE (> 3.0 A) - adjust protocol before docking compounds"

lines = HEAD + [
    "17 XCJ REDOCKING RMSD VALIDATION",
    "-" * 60,
    f"Reference     : {ref_path} (crystal pose, heavy atoms: {len(ref_xyz)})",
    f"Docked poses  : {poses_path} ({len(poses)} poses, export: {export_method})",
    "Method        : heavy-atom RMSD in the receptor frame (no superposition),",
    "                minimized over all symmetry-equivalent atom mappings.",
    "",
    f"{'POSE':<6}{'RMSD_A':>10}{'SYMMETRY_MAPS':>16}",
]
for n, r, m in rows:
    lines.append(f"{n:<6}{(f'{r:.3f}' if r is not None else 'NA'):>10}{(m or 0):>16}")

lines += [
    "",
    f"POSE_1_RMSD      {top:.3f} A" if top is not None else "POSE_1_RMSD      NA",
    f"LOWEST_RMSD      {best_row[1]:.3f} A (pose {best_row[0]})",
    "",
    "CRITERIA: < 2.0 A excellent | 2-3 A acceptable | > 3 A re-optimize",
    f"VERDICT: {verdict}",
]
(txt_dir / "17_XCJ_redock_RMSD_validation.txt").write_text("\n".join(lines) + "\n")
Path("XCJ_redock_verdict.txt").write_text(
    f"{top if top is not None else 'NA'}|{best_row[1]:.3f}|{best_row[0]}|{verdict}\n"
)
print("\n".join(lines[len(HEAD):]))
PY

# obrms cross-check (Open Babel, independent method)
if command -v obrms >/dev/null 2>&1; then
    {
        header_txt
        echo "17b OBRMS CROSS-CHECK (Open Babel, heavy-atom RMSD, symmetry-aware)"
        echo
        obrms "$XCJ_REF" XCJ_redock_poses.sdf 2>&1 | head -n 25
    } > "$TXT_DIR/17b_XCJ_redock_obrms_crosscheck.txt"
fi

stamp_sdf XCJ_redock_poses.sdf
stamp_mol2 XCJ_redock_best.mol2 "XCJ redocked best pose (Vina pose 1)"
stamp_mol2 XCJ_crystal.mol2 "XCJ crystal pose (7L00 chain ${SITE_CHAIN})"
stamp_pdb XCJ_redock_best.pdb "XCJ redocked best pose (Vina pose 1)"
stamp_pdb "$OUT" "XCJ redocking output (Vina)"

# ============================================================
# 2K. Step 2 summary
# ============================================================

BEST_SCORE=$(awk '/^[[:space:]]*1[[:space:]]+-?[0-9]/ {print $2; exit}' \
    "$TXT_DIR/16_XCJ_redock_scores.txt")
IFS='|' read -r RMSD1 RMSD_MIN RMSD_MIN_POSE VERDICT < XCJ_redock_verdict.txt

{
    header_txt
    echo "18 STEP 2 SUMMARY: RECEPTOR + REDOCKING VALIDATION"
    echo "------------------------------------------------------------"
    echo "Receptor           : ${RECEPTOR} (chains ${SITE_CHAIN}+${PARTNER_CHAIN} + FMN)"
    echo "Protein prep       : ${PROT_METHOD}"
    echo "Grid center        : ${CX}, ${CY}, ${CZ}"
    echo "Grid size          : ${SX} x ${SY} x ${SZ} A"
    echo "Exhaustiveness     : ${EXHAUSTIVENESS}   Modes: ${NUM_MODES}   Seed: ${SEED}"
    echo
    echo "XCJ best Vina score: ${BEST_SCORE:-NA} kcal/mol"
    echo "XCJ pose 1 RMSD    : ${RMSD1} A"
    echo "Lowest RMSD        : ${RMSD_MIN} A (pose ${RMSD_MIN_POSE})"
    echo "VERDICT            : ${VERDICT}"
    echo
    echo "Visual check (Discovery Studio): open ${TAG}_protein_FMN.pdb,"
    echo "XCJ_crystal.mol2 and XCJ_redock_best.mol2 together."
    echo
    echo "Vina scores are computational docking scores,"
    echo "not experimentally measured binding free energies."
} > "$TXT_DIR/18_step2_summary.txt"

rm -f XCJ_redock_verdict.txt "${TAG}_truncated_residues.txt"

echo
echo "============================================================"
echo "STEP 2 COMPLETED"
echo "============================================================"
cat "$TXT_DIR/18_step2_summary.txt"
