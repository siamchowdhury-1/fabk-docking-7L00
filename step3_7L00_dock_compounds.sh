#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 3: Dock pyrazole compounds into the validated 7L00 FabK
#         site (chains A+B + FMN), export complexes / MOL2,
#         analyse contacts, rank and graph.
#
# Usage:
#   ./step3_7L00_dock_compounds.sh                 # all *.mol in folder
#   ./step3_7L00_dock_compounds.sh 5sa23 17sa23    # selected compounds
#
# Optional: mic_values.txt  (one line per compound: name MIC_ug_per_mL)
#
# Focused-pocket run (small box on XCJ, separate folders):
#   RUN=focused BOX=16 ./step3_7L00_dock_compounds.sh
# Redraw graphs/ranking only (no docking):
#   GRAPHS_ONLY=1 ./step3_7L00_dock_compounds.sh
#   RUN=focused GRAPHS_ONLY=1 ./step3_7L00_dock_compounds.sh
# ============================================================

set -Eeuo pipefail

# ------------------------- settings -------------------------
PDB_ID="7L00"
TAG="7L00_AB"
RECEPTOR="${TAG}_receptor.pdbqt"
RECEPTOR_PDB="${TAG}_protein_FMN.pdb"
GRID="${TAG}_grid_box.txt"
PH="7.4"
EXHAUSTIVENESS="32"
NUM_MODES="20"
ENERGY_RANGE="5"
SEED="12345"
REF_DIR="txt_outputs"            # Step 2 reports live here
RUN="${RUN:-global}"             # global | focused | any label
BOX="${BOX:-}"                   # optional cubic box size override (A)
GRAPHS_ONLY="${GRAPHS_ONLY:-0}"

if [[ "$RUN" == "global" ]]; then
    TXT_DIR="txt_outputs"
    DIR_PREFIX="docking_"
    GRAPH_DIR="graphs"
else
    TXT_DIR="txt_outputs/${RUN}"
    DIR_PREFIX="docking_${RUN}_"
    GRAPH_DIR="graphs_${RUN}"
fi
# ------------------------------------------------------------

mkdir -p "$TXT_DIR"

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

stamp_pdb() {
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

stamp_mol2() {
    local f="$1" title="$2" tmp
    tmp="$(mktemp)"
    {
        for l in "${IDENTITY[@]}"; do echo "# $l"; done
        echo "# $title"
        grep -v '^# ' "$f" || true
    } > "$tmp"
    mv "$tmp" "$f"
}

stamp_sdf() {
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

# ------------------------- checks ---------------------------
for tool in python3 obabel vina mk_prepare_ligand.py; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not found. Run: conda activate docking"
        exit 1
    }
done

for f in "$RECEPTOR" "$RECEPTOR_PDB" "$GRID" \
         "$REF_DIR/08_${TAG}_active_site_residues.txt" \
         "$REF_DIR/18_step2_summary.txt"; do
    [[ -s "$f" ]] || { echo "ERROR: $f missing. Run Step 2 first."; exit 1; }
done

grep -q "VALIDATED\|ACCEPTABLE" "$REF_DIR/18_step2_summary.txt" || {
    echo "ERROR: Step 2 redocking did not pass validation."
    exit 1
}

read -r CX CY CZ SX SY SZ < <(awk '{printf "%s ", $2} END {print ""}' "$GRID")
if [[ -n "$BOX" ]]; then
    SX="$BOX"; SY="$BOX"; SZ="$BOX"
fi

if (( $# > 0 )); then
    COMPOUNDS=("$@")
else
    mapfile -t COMPOUNDS < <(ls *.mol 2>/dev/null | sed 's/\.mol$//' | sort -V)
fi

(( ${#COMPOUNDS[@]} > 0 )) || { echo "ERROR: no .mol files found."; exit 1; }

echo "============================================================"
echo "STEP 3 [${RUN}]: DOCKING ${#COMPOUNDS[@]} COMPOUND(S) INTO ${TAG}"
echo "Compounds: ${COMPOUNDS[*]}"
echo "Center: ${CX}, ${CY}, ${CZ}   Box: ${SX} x ${SY} x ${SZ} A"
echo "============================================================"

# ============================================================
# per-compound workflow
# ============================================================

dock_compound() {
    local cpd="$1"
    local dir="${DIR_PREFIX}${cpd}"
    local mol="${cpd}.mol"

    [[ -s "$mol" ]] || { echo "ERROR: $mol missing."; exit 1; }
    mkdir -p "$dir"

    echo
    echo "------------------------------------------------------------"
    echo "COMPOUND: ${cpd}"
    echo "------------------------------------------------------------"

    # ---------- 19 input check ----------
    {
        header_txt
        echo "19 INPUT CHECK: ${cpd}"
        echo "------------------------------------------------------------"
        ls -lh "$mol"
        echo
        echo "SMILES (as drawn):"
        obabel "$mol" -osmi 2>&1
        echo
        echo "MOL counts line:"
        sed -n '4p' "$mol"
        echo
        echo "Formula / MW (as drawn):"
        obabel "$mol" -otxt --append "formula MW" 2>/dev/null
        echo
        echo "Formal charge records:"
        grep '^M  CHG' "$mol" || echo "No explicit M  CHG record."
    } > "$TXT_DIR/19_${cpd}_input_check.txt"

    # ---------- 20 ligand preparation ----------
    local h_sdf="${dir}/${cpd}_H.sdf"
    local lig_qt="${dir}/${cpd}.pdbqt"
    local prep="$TXT_DIR/20_${cpd}_ligand_prepare.txt"

    rm -f "$h_sdf" "$lig_qt"
    {
        header_txt
        echo "20 LIGAND PREPARATION: ${cpd}"
        echo "------------------------------------------------------------"
        echo "obabel ${mol} -O ${h_sdf} -p ${PH} --gen3d --minimize --ff MMFF94 --steps 500"
        echo "mk_prepare_ligand.py -i ${h_sdf} -o ${lig_qt}"
        echo
    } > "$prep"

    obabel "$mol" -O "$h_sdf" -p "$PH" --gen3d \
        --minimize --ff MMFF94 --steps 500 >> "$prep" 2>&1
    [[ -s "$h_sdf" ]] || { echo "ERROR: 3D generation failed for $cpd"; exit 1; }

    mk_prepare_ligand.py -i "$h_sdf" -o "$lig_qt" >> "$prep" 2>&1
    [[ -s "$lig_qt" ]] || { echo "ERROR: PDBQT failed for $cpd (see $prep)"; exit 1; }

    python3 - "$h_sdf" "$SX" >> "$prep" <<'PY'
import math
import sys

path, box = sys.argv[1], float(sys.argv[2])
lines = open(path).read().splitlines()
n = int(lines[3][:3])
pts, heavy = [], 0
for l in lines[4:4 + n]:
    x, y, z = float(l[0:10]), float(l[10:20]), float(l[20:30])
    pts.append((x, y, z))
    if l[31:34].strip() != "H":
        heavy += 1
xs, ys, zs = zip(*pts)
span = [max(v) - min(v) for v in (xs, ys, zs)]
print()
print("3D DIMENSIONS (after protonation + MMFF94)")
print(f"ATOM_COUNT {n}   HEAVY_ATOMS {heavy}")
print(f"SIZE_X {span[0]:.3f} A   SIZE_Y {span[1]:.3f} A   SIZE_Z {span[2]:.3f} A")
print(f"MAXIMUM_DIMENSION {max(span):.3f} A   BOX {box:.1f} A")
print("FITS_BOX " + ("YES" if max(span) < box - 4 else "CHECK (ligand close to box size)"))
PY

    {
        echo
        echo "SMILES after protonation (pH ${PH}):"
        obabel "$h_sdf" -osmi 2>/dev/null
        echo "Net charge: $(awk '/^M  CHG/ {for (i=5; i<=NF; i+=2) q+=$i} /^M  END/ {print q+0; exit}' "$h_sdf")"
        echo "PDBQT atoms: $(grep -Ec '^(ATOM|HETATM)' "$lig_qt")"
        grep '^TORSDOF' "$lig_qt" || true
        echo
        echo "WARNINGS / ERRORS"
        grep -Ein 'error|warning|cannot|unknown|failed' "$prep" | grep -v 'WARNINGS / ERRORS' \
            || echo "None detected."
    } >> "$prep"

    stamp_sdf "$h_sdf"
    stamp_pdb "$lig_qt" "${cpd} ligand PDBQT (pH ${PH}, MMFF94, Meeko)"

    # ---------- 21 config + 22 vina ----------
    local config="${dir}/${cpd}_vina_config.txt"
    local out="${dir}/${cpd}_out.pdbqt"
    local vlog="$TXT_DIR/22_${cpd}_vina_run.txt"

    {
        for l in "${IDENTITY[@]}"; do echo "# $l"; done
        echo "# ${PDB_ID} FabK docking: ${cpd} (run: ${RUN}, box ${SX} A)"
        echo
        echo "receptor = ${RECEPTOR}"
        echo "ligand = ${lig_qt}"
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
    } > "$config"

    { header_txt; echo "21 VINA CONFIGURATION: ${cpd} (run ${RUN}, seed ${SEED})"; echo; cat "$config"; } \
        > "$TXT_DIR/21_${cpd}_vina_config.txt"

    { header_txt; echo "22 VINA RUN: ${cpd}"; echo; } > "$vlog"
    vina --config "$config" --out "$out" --seed "$SEED" 2>&1 | tee -a "$vlog"
    [[ -s "$out" ]] || { echo "ERROR: docking failed for $cpd"; exit 1; }

    {
        header_txt
        echo "23 DOCKING SCORES: ${cpd}"
        echo
        awk '/^mode[[:space:]]*\|/ {show=1} show {print}' "$vlog"
    } > "$TXT_DIR/23_${cpd}_scores.txt"

    # ---------- 24 export ----------
    local poses="${dir}/${cpd}_poses.sdf"
    local method="Meeko mk_export.py"
    rm -f "$poses"
    if command -v mk_export.py >/dev/null 2>&1; then
        mk_export.py "$out" -s "$poses" > /dev/null 2>&1 || true
    fi
    if [[ ! -s "$poses" ]]; then
        method="Open Babel (fallback, bond orders inferred)"
        obabel "$out" -O "$poses" 2>/dev/null
    fi

    obabel "$poses" -O "${dir}/${cpd}_all20.mol2" 2>/dev/null
    obabel "$poses" -l 1 -O "${dir}/${cpd}_best.mol2" 2>/dev/null
    obabel "$poses" -l 1 -O "${dir}/${cpd}_best_ligand.pdb" 2>/dev/null

    {
        header_txt
        echo "24 POSE EXPORT: ${cpd}"
        echo "Method: ${method}"
        echo
        ls -lh "$poses" "${dir}/${cpd}_all20.mol2" "${dir}/${cpd}_best.mol2" \
               "${dir}/${cpd}_best_ligand.pdb"
    } > "$TXT_DIR/24_${cpd}_pose_export.txt"

    # ---------- 25 complexes + contact analysis ----------
    python3 - "$cpd" "$poses" "$RECEPTOR_PDB" "$GRID" \
        "$REF_DIR/08_${TAG}_active_site_residues.txt" "$lig_qt" "$TXT_DIR" "$dir" <<'PY'
import math
import re
import sys
from datetime import datetime
from pathlib import Path

cpd, poses_path, rec_path, grid_path, site_path, lig_qt, txt_dir, out_dir = sys.argv[1:9]
txt_dir, out_dir = Path(txt_dir), Path(out_dir)

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
NOW = f"{datetime.now():%Y-%m-%d %H:%M:%S}"
HEAD = ["=" * 60, *IDENTITY, f"Date: {NOW}", "=" * 60, ""]


def read_sdf(path):
    records, block = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("$$$$"):
            records.append(block)
            block = []
        else:
            block.append(line)
    poses = []
    for rec in records:
        if len(rec) < 4:
            continue
        n = int(rec[3][:3])
        atoms = []
        for l in rec[4:4 + n]:
            atoms.append({
                "el": l[31:34].strip(),
                "xyz": (float(l[0:10]), float(l[10:20]), float(l[20:30])),
            })
        score = None
        for i, l in enumerate(rec):
            if re.match(r"^> *<(meeko|minimizedAffinity)", l) and i + 1 < len(rec):
                m = re.search(r'"energies":\s*\{[^}]*"free_energy":\s*(-?[\d.]+)', rec[i + 1])
                if m:
                    score = float(m.group(1))
                else:
                    try:
                        score = float(rec[i + 1].split()[0])
                    except ValueError:
                        pass
        poses.append({"atoms": atoms, "score": score})
    return poses


def read_receptor(path):
    atoms = []
    for l in Path(path).read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name = l[12:16].strip()
        el = l[76:78].strip() if len(l) >= 78 else ""
        if not el:
            el = re.sub(r"[^A-Za-z]", "", name)[:1]
        if el == "H":
            continue
        atoms.append({
            "line": l,
            "name": name,
            "res": l[17:20].strip(),
            "chain": l[21],
            "num": l[22:26].strip(),
            "el": el.upper(),
            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
        })
    return atoms


poses = read_sdf(poses_path)
rec = read_receptor(rec_path)
if not poses:
    raise SystemExit(f"ERROR: no poses read from {poses_path}")

# Vina scores from the docked PDBQT (authoritative)
out_pdbqt = Path(poses_path).with_name(f"{cpd}_out.pdbqt")
vina_scores = [
    float(l.split()[3]) for l in out_pdbqt.read_text().splitlines()
    if l.startswith("REMARK VINA RESULT:")
]
for i, p in enumerate(poses):
    if i < len(vina_scores):
        p["score"] = vina_scores[i]

grid = dict(l.split() for l in Path(grid_path).read_text().splitlines() if l.strip())
center = tuple(float(grid[k]) for k in ("CENTER_X", "CENTER_Y", "CENTER_Z"))

site_res = set()
for l in Path(site_path).read_text().splitlines():
    m = re.match(r"^([A-Za-z0-9])\s+([A-Z]{3})\s+(-?\d+\S*)\s+[\d.]+\s", l)
    if m:
        site_res.add((m.group(1), m.group(3), m.group(2)))

torsdof = "NA"
for l in Path(lig_qt).read_text().splitlines():
    if l.startswith("TORSDOF"):
        torsdof = l.split()[1]


def hetatm_lines(atoms, pose_no, start):
    out, counts = [], {}
    for k, a in enumerate(atoms):
        counts[a["el"]] = counts.get(a["el"], 0) + 1
        name = f"{a['el']}{counts[a['el']]}"[:4]
        x, y, z = a["xyz"]
        out.append(
            f"HETATM{(start + k) % 100000:5d} {name:<4} LIG Z{pose_no:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a['el']:>2}"
        )
    return out


def write_complex(path, pose_list, title):
    lines = [f"REMARK 999 {t}" for t in IDENTITY]
    lines += [f"REMARK 999 {title}", f"REMARK 999 Created: {NOW}"]
    lines += [a["line"] for a in rec] + ["TER"]
    serial = len(rec) + 1
    for n, p in pose_list:
        lines.append(f"REMARK 999 POSE {n} VINA_SCORE {p['score']}")
        lines += hetatm_lines(p["atoms"], n, serial)
        serial += len(p["atoms"])
        lines.append("TER")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


write_complex(
    out_dir / f"{cpd}_7L00_best_complex.pdb",
    [(1, poses[0])],
    f"7L00 FabK (A+B+FMN) + {cpd} best pose",
)
write_complex(
    out_dir / f"{cpd}_7L00_all_poses_complex.pdb",
    list(enumerate(poses, start=1)),
    f"7L00 FabK (A+B+FMN) + all {len(poses)} {cpd} poses (overlay, not simultaneous)",
)

# ---------- contact analysis for pose 1 ----------
best = poses[0]
lig = [a for a in best["atoms"] if a["el"] != "H"]
heavy = len(lig)
cent = tuple(sum(a["xyz"][i] for a in lig) / heavy for i in range(3))
cent_shift = math.dist(cent, center)

contacts, polar, hydrophobic, fmn_contacts = {}, [], set(), []
for r in rec:
    for a in lig:
        d = math.dist(r["xyz"], a["xyz"])
        if d > 4.0:
            continue
        key = (r["chain"], r["num"], r["res"])
        contacts[key] = min(contacts.get(key, 99.0), d)
        if a["el"] in ("N", "O") and r["el"] in ("N", "O") and d <= 3.5:
            polar.append((key, r["name"], a["el"], d))
        if a["el"] == "C" and r["el"] == "C" and r["res"] != "FMN":
            hydrophobic.add(key)
        if r["res"] == "FMN":
            fmn_contacts.append((r["chain"], r["name"], a["el"], d))

shared = sorted({k for k in contacts if k in site_res})
le = -best["score"] / heavy if best["score"] is not None else float("nan")

lines = HEAD + [
    f"25 BEST-POSE COMPLEX AND CONTACT ANALYSIS: {cpd}",
    "-" * 60,
    f"Best Vina score            : {best['score']} kcal/mol",
    f"Heavy atoms                : {heavy}",
    f"Ligand efficiency          : {le:.3f} kcal/mol per heavy atom",
    f"TORSDOF (rotatable bonds)  : {torsdof}",
    f"Pose-1 centroid to site    : {cent_shift:.2f} A from crystal XCJ center",
    f"Residues within 4.0 A      : {len(contacts)}",
    f"Shared with XCJ pocket     : {len(shared)} of {len(site_res)} (08 list, 5 A)",
    f"FMN contacts (<= 4.0 A)    : {'YES' if fmn_contacts else 'NO'} ({len(fmn_contacts)} atom pairs)",
    f"Polar contacts (N/O <=3.5A): {len(polar)} (possible H-bonds)",
    "",
    "RESIDUES WITHIN 4.0 A OF POSE 1",
    f"{'CHAIN':<7}{'RES':<6}{'NUM':<7}{'MIN_A':>7}  {'XCJ_SITE':<9}{'HYDROPHOBIC'}",
]
for key in sorted(contacts, key=lambda k: contacts[k]):
    ch, num, res = key
    lines.append(
        f"{ch:<7}{res:<6}{num:<7}{contacts[key]:>7.2f}  "
        f"{'yes' if key in site_res else '-':<9}{'yes' if key in hydrophobic else '-'}"
    )

lines += ["", "POLAR CONTACTS (ligand N/O to receptor N/O, <= 3.5 A)"]
lines += [
    f"  {k[0]}:{k[2]}{k[1]:<6} {rn:<5} ... ligand {le_} {d:.2f} A"
    for k, rn, le_, d in sorted(polar, key=lambda x: x[3])
] or ["  none"]

lines += ["", "FMN CONTACTS (<= 4.0 A)"]
lines += [
    f"  FMN {c}:{n:<5} ... ligand {e} {d:.2f} A"
    for c, n, e, d in sorted(fmn_contacts, key=lambda x: x[3])[:15]
] or ["  none"]

lines += [
    "",
    "ALL POSES",
    f"{'POSE':<6}{'SCORE':>8}{'CENTROID_TO_SITE_A':>21}",
]
for n, p in enumerate(poses, start=1):
    h = [a for a in p["atoms"] if a["el"] != "H"]
    c = tuple(sum(a["xyz"][i] for a in h) / len(h) for i in range(3))
    lines.append(f"{n:<6}{p['score']:>8}{math.dist(c, center):>21.2f}")

lines += [
    "",
    "FILES",
    f"  {out_dir}/{cpd}_7L00_best_complex.pdb",
    f"  {out_dir}/{cpd}_7L00_all_poses_complex.pdb",
    f"  {out_dir}/{cpd}_best.mol2  (Discovery Studio, with {rec_path})",
    "",
    "Geometric contacts only; use Discovery Studio for a full",
    "2D interaction diagram (H-bond angles, pi-stacking, salt bridges).",
]
(txt_dir / f"25_{cpd}_complex_contacts.txt").write_text("\n".join(lines) + "\n")

Path(out_dir / f"{cpd}_summary.tsv").write_text(
    "\t".join(str(v) for v in [
        cpd, best["score"], heavy, f"{le:.3f}", torsdof, f"{cent_shift:.2f}",
        len(contacts), len(shared), len(site_res),
        "YES" if fmn_contacts else "NO", len(polar),
    ]) + "\n"
)
print(f"{cpd}: best {best['score']} kcal/mol, {len(contacts)} contact residues, "
      f"FMN contact {'YES' if fmn_contacts else 'NO'}")
PY

    stamp_sdf "$poses"
    stamp_mol2 "${dir}/${cpd}_all20.mol2" "${cpd} all docked poses (7L00 FabK)"
    stamp_mol2 "${dir}/${cpd}_best.mol2" "${cpd} best docked pose (7L00 FabK)"
    stamp_pdb "${dir}/${cpd}_best_ligand.pdb" "${cpd} best docked pose (ligand only)"
    stamp_pdb "$out" "${cpd} Vina docking output (7L00 FabK)"
}

if [[ "$GRAPHS_ONLY" != "1" ]]; then
    for cpd in "${COMPOUNDS[@]}"; do
        dock_compound "$cpd"
    done
fi

# ============================================================
# 26 ranking + 27 graphs
# ============================================================

XCJ_SCORE=$(awk '/XCJ best Vina score/ {print $(NF-1)}' "$REF_DIR/18_step2_summary.txt")

python3 - "$TXT_DIR" "${XCJ_SCORE:-NA}" "$DIR_PREFIX" "$GRAPH_DIR" "$RUN" "$SX" \
    "${COMPOUNDS[@]}" <<'PY'
import math
import sys
from datetime import datetime
from pathlib import Path

txt_dir = Path(sys.argv[1])
xcj = sys.argv[2]
prefix, graph_dir, run, box = sys.argv[3:7]
cpds = sys.argv[7:]
label = "" if run == "global" else f" [{run}, {box} A box]"

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
NOW = f"{datetime.now():%Y-%m-%d %H:%M:%S}"
HEAD = ["=" * 60, *IDENTITY, f"Date: {NOW}", "=" * 60, ""]
FOOTER = "siam chowdhury | Computational and Medicinal Chemistry | " \
         "Dr. Alam's Research Team (www.alamresearch.org) | Arkansas State University"

rows = []
for c in cpds:
    f = Path(f"{prefix}{c}/{c}_summary.tsv")
    v = f.read_text().split("\t")
    rows.append({
        "cpd": v[0], "score": float(v[1]), "heavy": int(v[2]), "le": float(v[3]),
        "tors": v[4], "shift": float(v[5]), "contacts": int(v[6]),
        "shared": int(v[7]), "site": int(v[8]), "fmn": v[9], "polar": int(v[10].strip()),
    })

mic = {}
mic_file = Path("mic_values.txt")
if mic_file.is_file():
    for l in mic_file.read_text().splitlines():
        p = l.split()
        if len(p) >= 2 and not l.startswith("#"):
            try:
                mic[p[0]] = float(p[1])
            except ValueError:
                pass

rows.sort(key=lambda r: r["score"])

lines = HEAD + [
    f"26 DOCKING RANKING: 7L00 FabK (chains A+B + FMN){label}",
    "-" * 60,
    f"Reference: crystal inhibitor XCJ redocked = {xcj} kcal/mol (RMSD-validated)",
    "",
    f"{'RANK':<5}{'COMPOUND':<11}{'SCORE':>7}{'dXCJ':>7}{'LE':>7}{'HEAVY':>6}"
    f"{'TORS':>5}{'SHIFT_A':>8}{'CONT':>5}{'SITE':>7}{'FMN':>5}{'POLAR':>6}"
    + (f"{'MIC':>7}" if mic else ""),
]
for i, r in enumerate(rows, start=1):
    d = f"{r['score'] - float(xcj):+.2f}" if xcj != "NA" else "NA"
    lines.append(
        f"{i:<5}{r['cpd']:<11}{r['score']:>7.2f}{d:>7}{r['le']:>7.3f}{r['heavy']:>6}"
        f"{r['tors']:>5}{r['shift']:>8.2f}{r['contacts']:>5}"
        f"{str(r['shared']) + '/' + str(r['site']):>7}{r['fmn']:>5}{r['polar']:>6}"
        + (f"{mic.get(r['cpd'], float('nan')):>7.2f}" if mic else "")
    )

lines += [
    "",
    "COLUMNS",
    "  SCORE  best Vina score (kcal/mol, more negative = stronger predicted binding)",
    "  dXCJ   score minus XCJ reference (negative = better than XCJ)",
    "  LE     ligand efficiency, -score / heavy atoms",
    "  SHIFT  pose-1 centroid distance from crystal XCJ center",
    "  CONT   residues within 4 A;  SITE  shared with XCJ pocket residues",
    "  FMN    pose 1 contacts FMN;  POLAR  N/O pairs within 3.5 A",
]

pairs = [(r["score"], mic[r["cpd"]]) for r in rows if r["cpd"] in mic]
if len(pairs) >= 3:
    xs = [p[1] for p in pairs]
    ys = [p[0] for p in pairs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    r_val = sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")
    lines += [
        "",
        f"MIC vs SCORE (n={len(pairs)}): Pearson r = {r_val:.3f}",
        "Expected trend: lower MIC (more potent) <-> more negative score (r > 0).",
    ]
elif mic:
    lines += ["", "MIC correlation needs at least 3 compounds with MIC values."]

if run != "global":
    cmp_rows = []
    for r in rows:
        g = Path(f"docking_{r['cpd']}/{r['cpd']}_summary.tsv")
        if g.is_file():
            v = g.read_text().split("\t")
            cmp_rows.append((r, float(v[1]), float(v[5]), v[9]))
    if cmp_rows:
        lines += [
            "",
            f"COMPARISON WITH GLOBAL RUN (22 A box)",
            f"{'COMPOUND':<11}{'GLOBAL':>8}{'SHIFT':>7}{'FMN':>5}   "
            f"{run.upper():>8}{'SHIFT':>7}{'FMN':>5}{'PENALTY':>9}",
        ]
        for r, gs, gsh, gf in cmp_rows:
            lines.append(
                f"{r['cpd']:<11}{gs:>8.2f}{gsh:>7.2f}{gf:>5}   "
                f"{r['score']:>8.2f}{r['shift']:>7.2f}{r['fmn']:>5}"
                f"{r['score'] - gs:>+9.2f}"
            )
        lines += [
            "",
            "PENALTY = focused score - global score (kcal/mol).",
            "Small penalty (< ~1) with SHIFT < ~3 A and FMN contact: the compound",
            "fits the XCJ/FMN pocket. Large penalty or SHIFT still > 3 A: it does",
            "not fit the rigid crystal pocket well.",
        ]

lines += [
    "",
    "Vina scores are computational docking scores,",
    "not experimentally measured binding free energies.",
]
(txt_dir / "26_docking_ranking_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines[len(HEAD):]))

# ---------------- graphs ----------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
except ImportError:
    (txt_dir / "27_graphs_report.txt").write_text(
        "\n".join(HEAD + ["matplotlib not installed: conda install -c conda-forge matplotlib"]) + "\n"
    )
    raise SystemExit(0)


def mode_scores(c):
    modes, aff = [], []
    for l in Path(txt_dir / f"23_{c}_scores.txt").read_text().splitlines():
        p = l.split()
        try:
            m, a = int(p[0]), float(p[1])
        except (ValueError, IndexError):
            continue
        modes.append(m)
        aff.append(a)
    return modes, aff


def finish(fig, path):
    fig.text(0.5, 0.005, FOOTER, ha="center", va="bottom", fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


made = []
Path(graph_dir).mkdir(exist_ok=True)

for c in cpds:
    m, a = mode_scores(c)
    fig = plt.figure(figsize=(10, 6))
    plt.plot(m, a, marker="o", linewidth=2, markersize=6)
    for x, y in zip(m, a):
        plt.annotate(f"{y:.2f}", (x, y), xytext=(0, 8), textcoords="offset points",
                     ha="center", fontsize=8)
    plt.xlabel("Docking Mode")
    plt.ylabel("Binding Affinity (kcal/mol)")
    plt.title(f"Docking Scores: 7L00 FabK vs {c}{label}")
    plt.xticks(m)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle="--", alpha=0.4)
    out = f"{graph_dir}/Graph_7L00_FabK_{c}.png"
    finish(fig, out)
    made.append(out)

fig = plt.figure(figsize=(11, 6))
for c in cpds:
    m, a = mode_scores(c)
    plt.plot(m, a, marker="o", linewidth=2, label=c)
if xcj != "NA":
    plt.axhline(float(xcj), color="black", linestyle="--", linewidth=1,
                label=f"XCJ reference ({float(xcj):.2f})")
plt.xlabel("Docking Mode")
plt.ylabel("Binding Affinity (kcal/mol)")
plt.title(f"Docking Scores: 7L00 FabK vs Pyrazole Compounds{label}")
plt.xticks(range(1, 21))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
plt.gca().invert_yaxis()
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
out = f"{graph_dir}/multi_docking_7L00_FabK.png"
finish(fig, out)
made.append(out)

names = [r["cpd"] for r in rows]
scores = [r["score"] for r in rows]
if xcj != "NA":
    names.append("XCJ (ref)")
    scores.append(float(xcj))
fig = plt.figure(figsize=(max(6, 1.2 * len(names) + 3), 6))
colors = ["#4C72B0"] * len(rows) + (["#999999"] if xcj != "NA" else [])
bars = plt.bar(names, scores, color=colors)
for b, s in zip(bars, scores):
    plt.annotate(f"{s:.2f}", (b.get_x() + b.get_width() / 2, s), xytext=(0, 4),
                 textcoords="offset points", ha="center", va="bottom", fontsize=9)
plt.ylabel("Best Binding Affinity (kcal/mol)")
plt.title(f"Best Docking Score per Compound: 7L00 FabK{label}")
plt.ylim(0, min(scores) * 1.12)   # 0 at bottom, most negative at top
plt.grid(True, axis="y", linestyle="--", alpha=0.4)
out = f"{graph_dir}/best_scores_7L00_FabK.png"
finish(fig, out)
made.append(out)

if len(pairs) >= 2:
    fig = plt.figure(figsize=(7, 6))
    for r in rows:
        if r["cpd"] in mic:
            plt.scatter(mic[r["cpd"]], r["score"], s=60)
            plt.annotate(r["cpd"], (mic[r["cpd"]], r["score"]), xytext=(5, 5),
                         textcoords="offset points", fontsize=9)
    plt.xlabel("MIC (ug/mL)")
    plt.ylabel("Best Binding Affinity (kcal/mol)")
    plt.title("MIC vs Docking Score: 7L00 FabK")
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle="--", alpha=0.4)
    out = f"{graph_dir}/MIC_vs_score_7L00_FabK.png"
    finish(fig, out)
    made.append(out)

(txt_dir / "27_graphs_report.txt").write_text(
    "\n".join(HEAD + ["27 GRAPHS CREATED", "-" * 60] + made) + "\n"
)
print("\nGraphs:\n  " + "\n  ".join(made))
PY

echo
echo "============================================================"
echo "STEP 3 [${RUN}] COMPLETED"
echo "============================================================"
echo "Ranking : cat $TXT_DIR/26_docking_ranking_summary.txt"
echo "Contacts: cat $TXT_DIR/25_<compound>_complex_contacts.txt"
echo "Graphs  : ${GRAPH_DIR}/"
echo "Files   : ${DIR_PREFIX}<compound>/"
