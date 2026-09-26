#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 6: MM/GBSA VALIDATION CHECKS
#
# Step 4 ranked 20sa23 first, against both the Vina scores and
# the contact analysis. Every EEL term came out positive, which
# means electrostatics opposed binding in all three cases. That
# is the signature of anion-anion repulsion between the ligand
# carboxylate and Glu136. Three checks decide whether the
# ranking is real:
#
#   6A  XCJ control      MM/GBSA on the redocked XCJ pose, which
#                        we know is right to 0.82 A. If a correct
#                        pose scores badly, the protocol is at fault.
#   6B  His143 state     tleap made every histidine neutral (HIE).
#                        His143 is catalytic and may carry +1 (HIP),
#                        which would reverse the local electrostatics.
#   6C  Decomposition    per-residue energies: which residues help
#                        and which hurt, Glu136 above all.
#
# Usage:
#   ./step6_mmgbsa_checks.sh                     XCJ + all compounds, HIE and HIP
#   STATES="HIE" ./step6_mmgbsa_checks.sh        one protonation state only
#   MAXCYC=300 ./step6_mmgbsa_checks.sh          faster minimisation
#   ./step6_mmgbsa_checks.sh 5sa23               one compound (XCJ always included)
# ============================================================

set -Eeuo pipefail

# ------------------------- settings -------------------------
PDB_ID="7L00"
TAG="7L00_AB"
RUN="${RUN:-focused}"
CHAINS="${CHAINS:-A}"
CATALYTIC_HIS="${CATALYTIC_HIS:-143}"
STATES="${STATES:-HIE HIP}"    # HIE neutral (Nε-H), HID neutral (Nδ-H), HIP +1
FMN_CHARGE="${FMN_CHARGE:--2}"
MAXCYC="${MAXCYC:-500}"
NCYC="${NCYC:-200}"
IGB="${IGB:-5}"
SALT="${SALT:-0.15}"
DECOMP_CUT="${DECOMP_CUT:-6}"  # residues within this distance get decomposed
STEP4_DIR="mmgbsa_${RUN}"      # topologies and ligand parameters from Step 4
# ------------------------------------------------------------

POSE_PREFIX=$([[ "$RUN" == "global" ]] && echo "docking_" || echo "docking_${RUN}_")
TXT_DIR=$([[ "$RUN" == "global" ]] && echo "txt_outputs" || echo "txt_outputs/${RUN}")
WORK="mmgbsa_checks_${RUN}"

mkdir -p "$WORK" "$TXT_DIR"

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

hms() {
    local t="$1"
    if   (( t >= 3600 )); then printf '%dh %02dm %02ds' $((t/3600)) $((t%3600/60)) $((t%60))
    elif (( t >= 60 ));   then printf '%dm %02ds' $((t/60)) $((t%60))
    else                       printf '%ds' "$t"
    fi
}

begin() { STEP_LABEL="$1"; STEP_START=$SECONDS; printf '    %-36s' "${STEP_LABEL} ..."; }
done_() { local t=$(( SECONDS - STEP_START )); printf 'done in %s\n' "$(hms "$t")"; }
progress() {
    local pid="$1" start=$SECONDS
    while kill -0 "$pid" 2>/dev/null; do
        sleep 15
        printf '\r    %-36s%s elapsed' "${STEP_LABEL} ..." "$(hms $((SECONDS - start)))"
    done
    printf '\r'
}

stamp_pdb() {
    local f="$1" title="$2" tmp
    [[ -s "$f" ]] || return 0
    tmp="$(mktemp)"
    {
        for l in "${IDENTITY[@]}"; do echo "REMARK 999 $l"; done
        echo "REMARK 999 $title"
        grep -v '^REMARK 999 ' "$f" || true
    } > "$tmp"
    mv "$tmp" "$f"
}

# ------------------------- checks ---------------------------
for tool in antechamber parmchk2 tleap sander MMPBSA.py cpptraj obabel python3; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not found. Run: conda activate docking"; exit 1; }
done

for f in "$STEP4_DIR/protein_amber.pdb" "$STEP4_DIR/FMN.mol2" "$STEP4_DIR/FMN.frcmod"; do
    [[ -s "$f" ]] || { echo "ERROR: $f missing. Run Step 4 first."; exit 1; }
done
[[ -s "XCJ_redock_poses.sdf" ]] || {
    echo "ERROR: XCJ_redock_poses.sdf missing. Run Step 2 first."; exit 1; }

if (( $# > 0 )); then
    COMPOUNDS=("$@")
else
    mapfile -t COMPOUNDS < <(
        for d in ${POSE_PREFIX}*/; do
            [[ -d "$d" ]] && basename "$d" | sed "s/^${POSE_PREFIX}//"
        done | sort -V
    )
fi

LIGANDS=("XCJ" "${COMPOUNDS[@]}")

echo "============================================================"
echo "STEP 6: MM/GBSA VALIDATION CHECKS"
echo "Ligands  : ${LIGANDS[*]}"
echo "States   : ${STATES}  (His${CATALYTIC_HIS} of chain ${CHAINS:0:1})"
echo "Runs     : $(( ${#LIGANDS[@]} * $(wc -w <<< "$STATES") )) systems"
echo "============================================================"

# ============================================================
# 6-0. Locate His143 in the renumbered Amber PDB
# ============================================================

HIS_INDEX=$(python3 - "$STEP4_DIR" "$CHAINS" "$CATALYTIC_HIS" <<'PY'
import re
import sys
from pathlib import Path

work, chains, target = sys.argv[1], sys.argv[2][:1], sys.argv[3]

# pdb4amber renumbers from 1 and writes the old -> new map.
renum = Path(work) / "protein_amber_renum.txt"
hit = None
if renum.is_file():
    for line in renum.read_text(errors="replace").splitlines():
        f = line.split()
        if len(f) >= 5 and f[1] == chains and f[2] == target:
            hit = f[-1]
            break
        if len(f) >= 4 and f[1] == target and f[-2] in ("HIS", "HIE", "HID", "HIP"):
            hit = f[-1]

if hit is None:
    # Fall back to matching by residue name and order in the file itself.
    pdb = Path(work) / "protein_amber.pdb"
    seen = []
    for line in pdb.read_text(errors="replace").splitlines():
        if line.startswith("ATOM") and line[17:20].strip() in ("HIS", "HIE", "HID", "HIP"):
            num = line[22:26].strip()
            if num not in seen:
                seen.append(num)
    hit = ""

print(hit)
PY
)

if [[ -z "$HIS_INDEX" ]]; then
    echo "WARNING: could not map His${CATALYTIC_HIS} to the renumbered PDB."
    echo "         Only the default (HIE) state will be tested."
    STATES="HIE"
else
    echo "His${CATALYTIC_HIS} (chain ${CHAINS:0:1}) is residue ${HIS_INDEX} in protein_amber.pdb"
fi

# ============================================================
# 6-1. Receptor PDB for each protonation state
# ============================================================

for state in $STATES; do
    src="$STEP4_DIR/protein_amber.pdb"
    dst="$WORK/protein_${state}.pdb"

    if [[ "$state" == "HIE" || -z "$HIS_INDEX" ]]; then
        cp "$src" "$dst"      # tleap treats HIS as HIE by default
    else
        python3 - "$src" "$dst" "$HIS_INDEX" "$state" <<'PY'
import sys
from pathlib import Path

src, dst, index, state = sys.argv[1:5]
out, n = [], 0
for line in Path(src).read_text(errors="replace").splitlines():
    if (line.startswith(("ATOM", "HETATM"))
            and line[22:26].strip() == index
            and line[17:20].strip() in ("HIS", "HIE", "HID", "HIP")):
        line = line[:17] + f"{state:<3}" + line[20:]
        n += 1
    out.append(line)
Path(dst).write_text("\n".join(out) + "\n")
print(f"  {state}: {n} atoms renamed")
PY
    fi
done

# ============================================================
# 6-2. XCJ ligand parameters (redocked pose 1, the validated pose)
# ============================================================

XCJ_DIR="$WORK/XCJ_params"
mkdir -p "$XCJ_DIR"

if [[ ! -s "$XCJ_DIR/lig.mol2" ]]; then
    echo
    echo "--- XCJ control ligand ---"
    obabel XCJ_redock_poses.sdf -l 1 -O "$XCJ_DIR/lig_in.sdf" \
        > "$XCJ_DIR/obabel.log" 2>&1

    XCJ_NC=$(python3 - "$XCJ_DIR/lig_in.sdf" <<'PY'
import sys
try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = next(iter(Chem.SDMolSupplier(sys.argv[1], sanitize=False, removeHs=False)))
    print(sum(a.GetFormalCharge() for a in m.GetAtoms()))
except Exception:
    print(0)
PY
)
    echo "XCJ net charge: ${XCJ_NC}"

    begin "XCJ AM1-BCC charges"
    ( cd "$XCJ_DIR" && antechamber -i lig_in.sdf -fi sdf -o lig.mol2 -fo mol2 \
        -c bcc -nc "$XCJ_NC" -at gaff2 -rn LIG -s 2 -pf y ) \
        > "$XCJ_DIR/antechamber.log" 2>&1 &
    progress $!
    wait $! || true
    done_

    [[ -s "$XCJ_DIR/lig.mol2" ]] || {
        echo "ERROR: antechamber failed on XCJ."
        tail -n 15 "$XCJ_DIR/antechamber.log"
        exit 1
    }
    ( cd "$XCJ_DIR" && parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2 ) \
        > "$XCJ_DIR/parmchk2.log" 2>&1
fi

# ============================================================
# 6-3. One MM/GBSA system: ligand x protonation state
# ============================================================

run_system() {
    local lig="$1" state="$2"
    local dir="$WORK/${lig}_${state}"
    local src

    if [[ "$lig" == "XCJ" ]]; then
        src="$XCJ_DIR"
    else
        src="$STEP4_DIR/$lig"
    fi

    [[ -s "$src/lig.mol2" ]] || {
        echo "  ERROR: $src/lig.mol2 missing (Step 4 did not finish for $lig)."
        return 1
    }

    mkdir -p "$dir"
    cp "$src/lig.mol2" "$src/lig.frcmod" "$dir/"

    echo
    echo "  ---- ${lig} / His${CATALYTIC_HIS} = ${state} ----"

    # ---- topologies ----
    cat > "$dir/tleap.in" <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff2

loadamberparams ../../${STEP4_DIR}/FMN.frcmod
loadamberparams lig.frcmod

FMN = loadmol2 ../../${STEP4_DIR}/FMN.mol2
LIG = loadmol2 lig.mol2

prot = loadpdb ../protein_${state}.pdb

rec = combine { prot FMN }
com = combine { prot FMN LIG }

set default PBRadii mbondi2

saveamberparm LIG ligand.prmtop ligand.inpcrd
saveamberparm rec receptor.prmtop receptor.inpcrd
saveamberparm com complex.prmtop complex.inpcrd

quit
EOF

    begin "tleap"
    ( cd "$dir" && tleap -f tleap.in ) > "$dir/tleap.log" 2>&1 || true
    done_

    for f in complex.prmtop receptor.prmtop ligand.prmtop complex.inpcrd; do
        [[ -s "$dir/$f" ]] || {
            echo "    ERROR: tleap did not create $f"
            grep -Ei 'error|fatal|could not|unknown' "$dir/tleap.log" | head -n 10
            return 1
        }
    done

    # ---- minimisation ----
    cat > "$dir/min.in" <<EOF
Minimisation, GB implicit solvent, backbone restrained
 &cntrl
  imin = 1, maxcyc = ${MAXCYC}, ncyc = ${NCYC},
  ntb = 0, igb = ${IGB}, saltcon = ${SALT},
  cut = 999.0, ntpr = 100,
  ntr = 1, restraint_wt = 5.0,
  restraintmask = ':1-9999 & @CA,C,N,O',
 /
EOF

    begin "minimising ${MAXCYC} cycles"
    ( cd "$dir" && sander -O -i min.in -o min.out -p complex.prmtop \
        -c complex.inpcrd -ref complex.inpcrd -r complex_min.rst7 ) \
        > "$dir/sander.log" 2>&1 &
    progress $!
    wait $! || true
    done_

    [[ -s "$dir/complex_min.rst7" ]] || {
        echo "    ERROR: minimisation failed"
        tail -n 15 "$dir/min.out" 2>/dev/null
        return 1
    }

    ( cd "$dir" && cpptraj -p complex.prmtop <<'EOF'
trajin complex_min.rst7
trajout complex_min.nc netcdf
trajout complex_min.pdb pdb
go
quit
EOF
    ) > "$dir/cpptraj.log" 2>&1

    [[ -s "$dir/complex_min.nc" ]] || {
        echo "    ERROR: cpptraj failed"; tail -n 8 "$dir/cpptraj.log"; return 1; }

    # ---- residues near the ligand, for decomposition ----
    local res_list
    res_list=$(python3 - "$dir/complex_min.pdb" "$DECOMP_CUT" <<'PY'
import math
import sys
from pathlib import Path

pdb, cut = sys.argv[1], float(sys.argv[2])
prot, lig = [], []
for l in Path(pdb).read_text(errors="replace").splitlines():
    if not l.startswith(("ATOM", "HETATM")) or len(l) < 54:
        continue
    name = l[17:20].strip()
    xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    num = int(l[22:26])
    if name == "LIG":
        lig.append(xyz)
    else:
        prot.append((num, xyz))

near = set()
for num, p in prot:
    if num in near:
        continue
    for q in lig:
        if math.dist(p, q) <= cut:
            near.add(num)
            break

nums = sorted(near)
out, start, prev = [], None, None
for n in nums:
    if start is None:
        start = prev = n
    elif n == prev + 1:
        prev = n
    else:
        out.append(f"{start}-{prev}" if start != prev else f"{start}")
        start = prev = n
if start is not None:
    out.append(f"{start}-{prev}" if start != prev else f"{start}")
print(",".join(out))
PY
)

    # ---- MM/GBSA with per-residue decomposition ----
    cat > "$dir/mmgbsa.in" <<EOF
MM/GBSA with per-residue decomposition
&general
  startframe = 1, endframe = 1, interval = 1,
  verbose = 2, keep_files = 0,
/
&gb
  igb = ${IGB}, saltcon = ${SALT},
/
&decomp
  idecomp = 2, dec_verbose = 0,
  print_res = "${res_list}",
/
EOF

    begin "MM/GBSA + decomposition"
    ( cd "$dir" && MMPBSA.py -O -i mmgbsa.in -o mmgbsa.dat -eo frames.csv \
        -do decomp.dat -deo decomp_frames.csv \
        -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
        -y complex_min.nc ) > "$dir/mmpbsa.log" 2>&1 &
    progress $!
    wait $! || true
    done_

    if [[ ! -s "$dir/mmgbsa.dat" ]]; then
        echo "    decomposition failed, retrying without it"
        cat > "$dir/mmgbsa.in" <<EOF
MM/GBSA
&general
  startframe = 1, endframe = 1, interval = 1,
  verbose = 2, keep_files = 0,
/
&gb
  igb = ${IGB}, saltcon = ${SALT},
/
EOF
        begin "MM/GBSA (no decomposition)"
        ( cd "$dir" && MMPBSA.py -O -i mmgbsa.in -o mmgbsa.dat -eo frames.csv \
            -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
            -y complex_min.nc ) > "$dir/mmpbsa.log" 2>&1 &
        progress $!
        wait $! || true
        done_
    fi

    [[ -s "$dir/mmgbsa.dat" ]] || {
        echo "    ERROR: MMPBSA.py failed"; tail -n 20 "$dir/mmpbsa.log"; return 1; }

    stamp_pdb "$dir/complex_min.pdb" "${lig} complex, His${CATALYTIC_HIS}=${state}"
    echo "    OK: $(grep -E '^DELTA TOTAL' "$dir/mmgbsa.dat" | head -n1 | awk '{print $3}') kcal/mol"
}

RUN_START=$SECONDS
FAILED=()

for state in $STATES; do
    echo
    echo "============================================================"
    echo "PROTONATION STATE: His${CATALYTIC_HIS} = ${state}"
    echo "============================================================"
    for lig in "${LIGANDS[@]}"; do
        run_system "$lig" "$state" || FAILED+=("${lig}/${state}")
    done
done

TOTAL=$(( SECONDS - RUN_START ))

# ============================================================
# 6-4. Collect and interpret
# ============================================================

python3 - "$WORK" "$TXT_DIR" "$STEP4_DIR" "$POSE_PREFIX" "$CATALYTIC_HIS" \
    "$STATES" "${LIGANDS[@]}" <<'PY'
import math
import re
import sys
from datetime import datetime
from pathlib import Path

work, txt_dir, step4, prefix, his, states = sys.argv[1:7]
ligands = sys.argv[7:]
work, txt_dir = Path(work), Path(txt_dir)
states = states.split()

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
HEAD = ["=" * 76, *IDENTITY, f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}", "=" * 76, ""]
FOOTER = ("siam chowdhury | Computational and Medicinal Chemistry | "
          "Dr. Alam's Research Team (www.alamresearch.org) | Arkansas State University")

KEYS = ["VDWAALS", "EEL", "EGB", "ESURF", "DELTA TOTAL"]


def read_totals(path: Path):
    if not path.is_file():
        return None
    block = path.read_text(errors="replace").split(
        "Differences (Complex - Receptor - Ligand)")[-1]
    out = {}
    for k in KEYS:
        m = re.search(rf"^{re.escape(k)}\s+(-?\d+\.\d+)", block, re.M)
        if m:
            out[k] = float(m.group(1))
    return out if "DELTA TOTAL" in out else None


def read_decomp(path: Path):
    """Per-residue TOTAL contributions from the decomposition output."""
    if not path.is_file():
        return {}
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        # e.g.  GLU 136 | -1.234 | ... | TOTAL
        m = re.match(
            r"^\s*([A-Z]{2,4})\s+(\d+)\s*\|?.*?(-?\d+\.\d+)\s*$", line)
        if not m:
            continue
        f = [x for x in re.split(r"[,|]", line) if x.strip()]
        nums = [float(x) for x in re.findall(r"-?\d+\.\d+", line)]
        if len(nums) >= 2:
            out[f"{m.group(1)}{m.group(2)}"] = nums[-2] if len(nums) > 5 else nums[-1]
    return out


results = {}
for st in states:
    for lig in ligands:
        t = read_totals(work / f"{lig}_{st}" / "mmgbsa.dat")
        if t:
            results[(lig, st)] = {
                "totals": t,
                "decomp": read_decomp(work / f"{lig}_{st}" / "decomp.dat"),
            }

if not results:
    raise SystemExit("ERROR: no MM/GBSA results could be read from the checks.")

# Step 4 values for comparison (HIE, no decomposition)
step4_totals = {}
for lig in ligands:
    t = read_totals(Path(step4) / lig / "mmgbsa.dat")
    if t:
        step4_totals[lig] = t

vina = {}
for lig in ligands:
    f = Path(f"{prefix}{lig}/{lig}_summary.tsv")
    if f.is_file():
        vina[lig] = float(f.read_text().split("\t")[1])

lines = HEAD + [
    "39 MM/GBSA VALIDATION CHECKS",
    "-" * 76,
    "6A  XCJ control      the redocked crystal inhibitor, RMSD 0.82 A from",
    "                     its crystal pose. Its energy tests the protocol,",
    "                     because this pose is known to be right.",
    f"6B  His{his} state    HIE is neutral (proton on N-epsilon), HID neutral",
    "                     (N-delta), HIP protonated and positively charged.",
    "6C  Decomposition    per-residue contributions to the binding energy.",
    "",
    "=" * 76,
    "6A + 6B  TOTALS BY PROTONATION STATE",
    "=" * 76,
    f"{'LIGAND':<10}{'STATE':<7}{'MMGBSA':>9}{'VDW':>9}{'EEL':>10}{'EGB':>10}"
    f"{'ESURF':>8}{'VINA':>8}",
]

for st in states:
    for lig in ligands:
        r = results.get((lig, st))
        if not r:
            lines.append(f"{lig:<10}{st:<7}{'failed':>9}")
            continue
        t = r["totals"]
        v = vina.get(lig)
        lines.append(
            f"{lig:<10}{st:<7}{t['DELTA TOTAL']:>9.2f}{t.get('VDWAALS', 0):>9.2f}"
            f"{t.get('EEL', 0):>10.2f}{t.get('EGB', 0):>10.2f}{t.get('ESURF', 0):>8.2f}"
            f"{(v if v is not None else float('nan')):>8.2f}"
        )
    lines.append("")

# --- ranking per state ---
lines += ["=" * 76, "RANKING BY STATE (XCJ shown for reference)", "=" * 76]
for st in states:
    rows = sorted(
        [(l, results[(l, st)]["totals"]["DELTA TOTAL"])
         for l in ligands if (l, st) in results],
        key=lambda x: x[1],
    )
    order = " > ".join(f"{l} ({v:.1f})" for l, v in rows)
    lines.append(f"  His{his}={st:<4} {order}")

if vina:
    vrows = sorted(vina.items(), key=lambda x: x[1])
    lines.append("  Vina      " + " > ".join(f"{l} ({v:.1f})" for l, v in vrows))
lines.append("")

# --- XCJ control verdict ---
lines += ["=" * 76, "6A  XCJ CONTROL: IS THE PROTOCOL SOUND?", "=" * 76]
for st in states:
    r = results.get(("XCJ", st))
    if not r:
        lines.append(f"  His{his}={st}: XCJ run failed.")
        continue
    xt = r["totals"]["DELTA TOTAL"]
    comp = [results[(l, st)]["totals"]["DELTA TOTAL"]
            for l in ligands if l != "XCJ" and (l, st) in results]
    lines.append(f"  His{his}={st}: XCJ = {xt:.2f} kcal/mol, EEL = {r['totals'].get('EEL', 0):.2f}")
    if comp:
        better = sum(1 for c in comp if c < xt)
        lines.append(
            f"    {better} of {len(comp)} compounds score better than the known-correct pose."
        )
    if r["totals"].get("EEL", 0) > 0:
        lines.append(
            "    WARNING: electrostatics oppose binding even for the crystal"
            " inhibitor,"
        )
        lines.append(
            "    so the positive EEL values are a property of this setup, not of"
        )
        lines.append("    your compounds.")
    else:
        lines.append(
            "    EEL is favourable for the crystal inhibitor, so positive EEL"
        )
        lines.append("    for a compound is a real signal about that compound.")
lines.append("")

# --- state sensitivity ---
if len(states) > 1:
    lines += ["=" * 76, f"6B  DOES His{his} PROTONATION CHANGE THE ANSWER?", "=" * 76,
              f"{'LIGAND':<10}" + "".join(f"{s:>12}" for s in states) + f"{'SPREAD':>10}"]
    for lig in ligands:
        vals = [results[(lig, s)]["totals"]["DELTA TOTAL"] if (lig, s) in results
                else None for s in states]
        if any(v is None for v in vals):
            continue
        lines.append(
            f"{lig:<10}" + "".join(f"{v:>12.2f}" for v in vals)
            + f"{max(vals) - min(vals):>10.2f}"
        )

    orders = {}
    for st in states:
        rows = sorted(
            [(l, results[(l, st)]["totals"]["DELTA TOTAL"])
             for l in ligands if l != "XCJ" and (l, st) in results],
            key=lambda x: x[1],
        )
        orders[st] = [l for l, _ in rows]
    distinct = {tuple(v) for v in orders.values()}
    lines += [
        "",
        f"  Compound order is {'THE SAME' if len(distinct) == 1 else 'DIFFERENT'}"
        f" across protonation states.",
    ]
    if len(distinct) > 1:
        lines += [
            "  The ranking depends on an assumption about one proton, so it",
            "  cannot be reported as a result until that state is established.",
        ]
    lines.append("")

# --- decomposition ---
have_decomp = any(r["decomp"] for r in results.values())
if have_decomp:
    lines += ["=" * 76, "6C  PER-RESIDUE CONTRIBUTIONS (kcal/mol)", "=" * 76,
              "Negative helps binding, positive opposes it.", ""]
    for st in states:
        for lig in ligands:
            r = results.get((lig, st))
            if not r or not r["decomp"]:
                continue
            d = sorted(r["decomp"].items(), key=lambda x: x[1])
            n = min(6, len(d) // 2) or len(d)
            lines.append(f"  {lig} / His{his}={st}   ({len(d)} residues decomposed)")
            if len(d) <= 12:
                lines.append(f"    {'RESIDUE':<10}{'kcal/mol':>10}")
                for name, val in d:
                    lines.append(f"    {name:<10}{val:>10.2f}")
            else:
                lines.append(f"    {'most favourable':<26}{'most unfavourable'}")
                for i in range(n):
                    left = f"{d[i][0]:<10}{d[i][1]:>9.2f}"
                    j = len(d) - 1 - i
                    right = f"{d[j][0]:<10}{d[j][1]:>9.2f}"
                    lines.append(f"    {left:<26}{right}")
            lines.append("")

    # Glu136 specifically
    glu = [(lig, st, r["decomp"][k])
           for (lig, st), r in results.items()
           for k in r["decomp"] if k.startswith("GLU") and "136" in k]
    if glu:
        lines += ["  GLU136 CONTRIBUTION (the suspected repulsion)",
                  f"    {'LIGAND':<10}{'STATE':<7}{'kcal/mol':>10}"]
        for lig, st, v in sorted(glu, key=lambda x: (x[1], x[0])):
            lines.append(f"    {lig:<10}{st:<7}{v:>10.2f}")
        lines.append("")
else:
    lines += ["=" * 76, "6C  PER-RESIDUE CONTRIBUTIONS", "=" * 76,
              "  Decomposition did not run. See the mmpbsa.log files.", ""]

lines += [
    "=" * 76,
    "HOW TO READ THIS",
    "=" * 76,
    "Single-structure MM/GBSA with charged ligands is noisy. EEL and EGB are",
    "large numbers that nearly cancel, so a small error in either moves the",
    "total a long way. Treat differences under about 5 kcal/mol as undecided.",
    "",
    "No entropy term is included. These are relative ranking energies, not",
    "experimental binding free energies.",
    "",
]

(txt_dir / "39_mmgbsa_validation_checks.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines[len(HEAD):]))

# ---- graph: state sensitivity ----
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit(0)

graph_dir = Path("graphs" if prefix == "docking_" else "graphs_focused")
graph_dir.mkdir(exist_ok=True)

present = [l for l in ligands if any((l, s) in results for s in states)]
if present and len(states) >= 1:
    width = 0.8 / len(states)
    fig = plt.figure(figsize=(max(7, 1.7 * len(present) + 4), 6))
    for j, st in enumerate(states):
        vals = [results[(l, st)]["totals"]["DELTA TOTAL"] if (l, st) in results else 0
                for l in present]
        xs = [i + j * width - 0.4 for i in range(len(present))]
        bars = plt.bar(xs, vals, width=width, label=f"His{his}={st}")
        for b, v in zip(bars, vals):
            if v:
                plt.annotate(f"{v:.0f}", (b.get_x() + b.get_width() / 2, v),
                             xytext=(0, -12 if v < 0 else 4),
                             textcoords="offset points", ha="center", fontsize=8)
    plt.xticks(range(len(present)), present)
    plt.ylabel("MM/GBSA binding energy (kcal/mol)")
    plt.title(f"MM/GBSA by His{his} protonation state (XCJ = validated control)")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.legend()
    plt.grid(True, axis="y", linestyle="--", alpha=0.4)
    fig.text(0.5, 0.005, FOOTER, ha="center", va="bottom", fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(graph_dir / "mmgbsa_protonation_check.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nGraph: {graph_dir}/mmgbsa_protonation_check.png")
PY

{
    header_txt
    echo "40 STEP 6 RUN TIMINGS"
    echo "------------------------------------------------------------"
    echo "Settings: MAXCYC=${MAXCYC}, igb=${IGB}, states: ${STATES}"
    echo "Ligands : ${LIGANDS[*]}"
    echo "Total   : $(hms "$TOTAL")"
    if (( ${#FAILED[@]} )); then
        echo
        echo "Failed systems: ${FAILED[*]}"
    fi
} > "$TXT_DIR/40_step6_timings.txt"

echo
echo "============================================================"
echo "STEP 6 COMPLETED in $(hms "$TOTAL")"
echo "============================================================"
(( ${#FAILED[@]} )) && echo "FAILED: ${FAILED[*]}"
echo "Read: cat $TXT_DIR/39_mmgbsa_validation_checks.txt"
