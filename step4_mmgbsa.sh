#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 4: MM/GBSA rescoring of docked poses (AmberTools).
#
#   ligand + FMN parameters (GAFF2 / AM1-BCC)
#        -> tleap complex, receptor, ligand topologies (ff14SB)
#        -> sander minimisation (igb=5)
#        -> MMPBSA.py single-structure MM/GBSA
#        -> ranking vs Vina and anchor contacts
#
# Usage:
#   ./step4_mmgbsa.sh                          # focused run, all compounds
#   RUN=global ./step4_mmgbsa.sh               # global run poses
#   ./step4_mmgbsa.sh 5sa23 17sa23             # selected compounds
#   CHAINS=AB ./step4_mmgbsa.sh                # keep both chains (slower)
# ============================================================

set -Eeuo pipefail

# ------------------------- settings -------------------------
PDB_ID="7L00"
TAG="7L00_AB"
RUN="${RUN:-focused}"
CHAINS="${CHAINS:-A}"          # A = pocket chain only; AB = both
FMN_CHARGE="${FMN_CHARGE:--2}" # FMN phosphate, dianionic at pH 7.4
MAXCYC="${MAXCYC:-2000}"       # minimisation steps
NCYC="${NCYC:-500}"            # steepest-descent steps before conjugate gradient
IGB="${IGB:-5}"                # GB model for minimisation and MM/GBSA
SALT="${SALT:-0.15}"           # salt concentration, mol/L
REF_DIR="txt_outputs"
# ------------------------------------------------------------

if [[ "$RUN" == "global" ]]; then
    POSE_PREFIX="docking_"
    TXT_DIR="txt_outputs"
else
    POSE_PREFIX="docking_${RUN}_"
    TXT_DIR="txt_outputs/${RUN}"
fi

WORK="mmgbsa_${RUN}"
mkdir -p "$WORK" "$TXT_DIR"

IDENTITY=(
    "siam chowdhury"
    "Computational and Medicinal Chemistry"
    "[Dr. Alam's Research Team] www.alamresearch.org"
    "Arkansas State University"
)

# ---- timing ----
STEP_START=""

hms() {   # seconds -> 1h 02m 03s
    local t="$1"
    if   (( t >= 3600 )); then printf '%dh %02dm %02ds' $((t/3600)) $((t%3600/60)) $((t%60))
    elif (( t >= 60 ));   then printf '%dm %02ds' $((t/60)) $((t%60))
    else                       printf '%ds' "$t"
    fi
}

begin() {   # begin "label"
    STEP_LABEL="$1"
    STEP_START=$SECONDS
    printf '  %-34s' "${STEP_LABEL} ..."
}

done_() {   # print elapsed for the step just begun
    local t=$(( SECONDS - STEP_START ))
    printf 'done in %s\n' "$(hms "$t")"
    LAST_STEP_SECONDS=$t
}

progress() {   # background ticker so long steps show they are alive
    local pid="$1" start=$SECONDS
    while kill -0 "$pid" 2>/dev/null; do
        sleep 15
        printf '\r  %-34s%s elapsed' "${STEP_LABEL} ..." "$(hms $((SECONDS - start)))"
    done
    printf '\r'
}

header_txt() {
    echo "============================================================"
    printf '%s\n' "${IDENTITY[@]}"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"
    echo
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

# Net formal charge of a structure file. RDKit first (it reads the
# bond orders and formal charges directly), then M CHG records, then
# a SMILES count. A wrong charge here makes sqm fail with an odd
# electron count, so it is worth getting from more than one source.
net_charge() {
    python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path

path = sys.argv[1]
q = None

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = None
    if path.endswith((".sdf", ".mol")):
        supplier = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
    elif path.endswith(".mol2"):
        mol = Chem.MolFromMol2File(path, sanitize=False, removeHs=False)
    if mol is not None:
        q = sum(a.GetFormalCharge() for a in mol.GetAtoms())
except Exception:
    pass

if q is None:
    total, seen = 0, False
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith("M  CHG"):
            f = line.split()
            total += sum(int(f[i]) for i in range(5, len(f), 2))
            seen = True
        if line.startswith("M  END"):
            break
    if seen:
        q = total

if q is None:
    try:
        import subprocess
        smi = subprocess.run(
            ["obabel", path, "-osmi"], capture_output=True, text=True
        ).stdout
        q = len(re.findall(r"\[[^]]*-\]", smi)) * -1 + len(re.findall(r"\[[^]]*\+\]", smi))
    except Exception:
        q = 0

print(int(q))
PY
}

# ------------------------- checks ---------------------------
for tool in antechamber parmchk2 tleap sander MMPBSA.py cpptraj pdb4amber obabel python3; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not found. Run: conda activate docking"
        exit 1
    }
done

[[ -s "${TAG}_protein.pdb" ]] || { echo "ERROR: ${TAG}_protein.pdb missing (Step 2)."; exit 1; }
[[ -s "${PDB_ID}_FMN.pdb" ]]  || { echo "ERROR: ${PDB_ID}_FMN.pdb missing (Step 1)."; exit 1; }

if (( $# > 0 )); then
    COMPOUNDS=("$@")
else
    mapfile -t COMPOUNDS < <(
        for d in ${POSE_PREFIX}*/; do
            [[ -d "$d" ]] && basename "$d" | sed "s/^${POSE_PREFIX}//"
        done | sort -V
    )
fi
(( ${#COMPOUNDS[@]} > 0 )) || { echo "ERROR: no ${POSE_PREFIX}* folders found."; exit 1; }

echo "============================================================"
echo "STEP 4: MM/GBSA [${RUN} poses, chain(s) ${CHAINS}]"
echo "Compounds: ${COMPOUNDS[*]}"
echo "============================================================"

# ============================================================
# 4A. Receptor: protein chain(s) + FMN, cleaned for tleap
# ============================================================

echo
echo "--- receptor preparation ---"

python3 - "$TAG" "$PDB_ID" "$CHAINS" "$WORK" <<'PY'
import sys
from pathlib import Path

tag, pdb_id, chains, work = sys.argv[1:5]
work = Path(work)
keep = set(chains)

prot = [
    l for l in Path(f"{tag}_protein.pdb").read_text().splitlines()
    if l.startswith(("ATOM", "HETATM")) and l[21] in keep
]
fmn = [
    l for l in Path(f"{pdb_id}_FMN.pdb").read_text().splitlines()
    if l.startswith(("ATOM", "HETATM")) and l[21] in keep
]
if not prot:
    raise SystemExit(f"ERROR: no protein atoms for chain(s) {chains}")

# TER at chain changes and at numbering gaps, so tleap does not
# bond across a break in the crystal structure.
out, prev_ch, prev_num, gaps = [], None, None, []
for l in prot:
    ch, num = l[21], int(l[22:26])
    if prev_ch is not None and (ch != prev_ch or num - prev_num > 1):
        out.append("TER")
        if ch == prev_ch:
            gaps.append(f"{ch}:{prev_num}-{num}")
    out.append(l)
    prev_ch, prev_num = ch, num
out.append("TER")

Path(f"{work}/receptor_protein.pdb").write_text("\n".join(out) + "\nEND\n")
Path(f"{work}/FMN.pdb").write_text("\n".join(fmn) + "\nEND\n")
Path(f"{work}/receptor_notes.txt").write_text(
    f"chains {chains}\nprotein atoms {len(prot)}\nFMN atoms {len(fmn)}\n"
    f"chain breaks (TER inserted): {', '.join(gaps) if gaps else 'none'}\n"
)
print(f"protein atoms {len(prot)}, FMN atoms {len(fmn)}")
print(f"chain breaks: {', '.join(gaps) if gaps else 'none'}")
PY

pdb4amber -i "$WORK/receptor_protein.pdb" -o "$WORK/protein_amber.pdb" \
    --nohyd --dry > "$WORK/pdb4amber.log" 2>&1 || true
[[ -s "$WORK/protein_amber.pdb" ]] || {
    echo "ERROR: pdb4amber failed. See $WORK/pdb4amber.log"
    exit 1
}

# ============================================================
# 4B. FMN parameters (once, reused by every compound)
# ============================================================

FMN_SOURCE="crystal SDF from RCSB (CCD bond orders)"

if [[ ! -s "$WORK/FMN.mol2" ]]; then
    echo
    echo "--- FMN parameters (AM1-BCC, charge ${FMN_CHARGE}) ---"
    echo "This takes a few minutes."

    # Bond orders must come from a real chemical record. A PDB file has
    # none, and antechamber's bondtype cannot guess them around phosphate.
    FMN_SDF=""
    for cand in "FMN_${CHAINS:0:1}401_crystal.sdf" FMN_*_crystal.sdf; do
        if [[ -s "$cand" ]]; then FMN_SDF="$cand"; break; fi
    done

    if [[ -z "$FMN_SDF" ]]; then
        echo "FMN crystal SDF not found; downloading from RCSB..."
        url="https://models.rcsb.org/v1/${PDB_ID,,}/ligand?auth_asym_id=${CHAINS:0:1}&auth_seq_id=401&encoding=sdf"
        if command -v wget >/dev/null 2>&1; then
            wget -q -O "$WORK/FMN_rcsb.sdf" "$url" || true
        else
            curl -fsSL -o "$WORK/FMN_rcsb.sdf" "$url" || true
        fi
        [[ -s "$WORK/FMN_rcsb.sdf" ]] && FMN_SDF="$WORK/FMN_rcsb.sdf"
    fi

    if [[ -n "$FMN_SDF" ]]; then
        obabel "$FMN_SDF" -O "$WORK/FMN_h.mol2" -h -p 7.4 \
            > "$WORK/FMN_obabel.log" 2>&1
    else
        FMN_SOURCE="Open Babel perception from crystal PDB (fallback)"
        obabel "$WORK/FMN.pdb" -O "$WORK/FMN_h.mol2" -h -p 7.4 \
            > "$WORK/FMN_obabel.log" 2>&1
    fi

    [[ -s "$WORK/FMN_h.mol2" ]] || {
        echo "ERROR: could not build an FMN MOL2. See $WORK/FMN_obabel.log"
        exit 1
    }

    FMN_NC=$(net_charge "$WORK/FMN_h.mol2")
    echo "FMN source: ${FMN_SOURCE}"
    echo "FMN net charge from structure: ${FMN_NC} (setting: ${FMN_CHARGE})"

    if [[ "$FMN_NC" != "$FMN_CHARGE" ]]; then
        echo "NOTE: using the charge found in the structure (${FMN_NC})."
        FMN_CHARGE="$FMN_NC"
    fi

    # -j 1 = assign atom types only, keep the bond orders from the MOL2.
    begin "FMN AM1-BCC charges"
    ( cd "$WORK" && antechamber -i FMN_h.mol2 -fi mol2 -o FMN.mol2 -fo mol2 \
        -c bcc -nc "$FMN_CHARGE" -at gaff2 -rn FMN -j 1 -s 2 -pf y ) \
        > "$WORK/FMN_antechamber.log" 2>&1 &
    progress $!
    wait $! || true
    done_

    [[ -s "$WORK/FMN.mol2" ]] || {
        echo "ERROR: antechamber failed on FMN. See $WORK/FMN_antechamber.log"
        echo
        tail -n 20 "$WORK/FMN_antechamber.log"
        echo
        echo "If sqm did not converge, try: FMN_CHARGE=-1 ./step4_mmgbsa.sh"
        exit 1
    }
    ( cd "$WORK" && parmchk2 -i FMN.mol2 -f mol2 -o FMN.frcmod -s gaff2 ) \
        > "$WORK/FMN_parmchk2.log" 2>&1
fi

{
    header_txt
    echo "35 RECEPTOR AND FMN PARAMETERISATION"
    echo "------------------------------------------------------------"
    cat "$WORK/receptor_notes.txt"
    echo
    echo "Protein force field : ff14SB"
    echo "FMN / ligand        : GAFF2 with AM1-BCC charges"
    echo "FMN net charge      : ${FMN_CHARGE}"
    echo "FMN bond orders     : ${FMN_SOURCE:-reused from a previous run}"
    echo
    echo "pdb4amber changes (renamed / removed):"
    grep -Ei 'renam|remov|miss|gap|sum' "$WORK/pdb4amber.log" | head -n 30 \
        || echo "none reported"
    echo
    echo "FMN parameters missing from GAFF2 and estimated by parmchk2:"
    grep -c "ATTN" "$WORK/FMN.frcmod" 2>/dev/null || echo 0
} > "$TXT_DIR/35_mmgbsa_receptor_setup.txt"
cat "$TXT_DIR/35_mmgbsa_receptor_setup.txt"

# ============================================================
# 4C. Per-compound MM/GBSA
# ============================================================

run_compound() {
    local cpd="$1"
    local dir="$WORK/$cpd"
    local poses="${POSE_PREFIX}${cpd}/${cpd}_poses.sdf"

    [[ -s "$poses" ]] || { echo "ERROR: $poses missing."; return 1; }
    mkdir -p "$dir"

    echo
    echo "------------------------------------------------------------"
    echo "MM/GBSA: ${cpd}"
    echo "------------------------------------------------------------"

    # ---- best pose as SDF: bond orders AND formal charges survive ----
    obabel "$poses" -l 1 -O "$dir/lig_in.sdf" > "$dir/obabel.log" 2>&1
    obabel "$poses" -l 1 -O "$dir/lig_in.mol2" >> "$dir/obabel.log" 2>&1
    [[ -s "$dir/lig_in.sdf" ]] || { echo "ERROR: could not extract pose 1."; return 1; }

    # The SDF keeps formal charges, so it is the primary source. A MOL2
    # conversion can silently drop them. Where the two disagree, the
    # electron count decides: a closed-shell molecule must have an even
    # number, and sqm refuses to run when it does not.
    local nc nc_sdf nc_mol2 zsum
    nc_sdf=$(net_charge "$dir/lig_in.sdf")
    nc_mol2=$(net_charge "$dir/lig_in.mol2")
    zsum=$(python3 - "$dir/lig_in.sdf" <<'PY'
import sys
from pathlib import Path
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
     "CL": 17, "BR": 35, "I": 53}
lines = Path(sys.argv[1]).read_text(errors="replace").splitlines()
n = int(lines[3][:3])
print(sum(Z.get(l.split()[3].upper(), 0) for l in lines[4:4 + n]))
PY
)

    nc=""
    for cand in "$nc_sdf" "$nc_mol2"; do
        if (( (zsum - cand) % 2 == 0 )); then nc="$cand"; break; fi
    done

    if [[ -z "$nc" ]]; then
        echo "ERROR: ${cpd}: no candidate charge gives an even electron count."
        echo "  electrons at neutral: ${zsum}"
        echo "  charge from SDF ${nc_sdf}, from MOL2 ${nc_mol2}"
        echo "Check the structure: obabel ${dir}/lig_in.sdf -osmi"
        return 1
    fi

    if [[ "$nc_sdf" != "$nc_mol2" ]]; then
        echo "NOTE: charge from SDF ${nc_sdf}, from MOL2 ${nc_mol2};"
        echo "      using ${nc} (gives an even electron count)."
    fi
    echo "net charge: ${nc}   electrons: $((zsum - nc))"

    # ---- GAFF2 / AM1-BCC parameters ----
    local t_start=$SECONDS

    if [[ ! -s "$dir/lig.mol2" ]]; then
        begin "AM1-BCC charges"
        ( cd "$dir" && antechamber -i lig_in.sdf -fi sdf -o lig.mol2 -fo mol2 \
            -c bcc -nc "$nc" -at gaff2 -rn LIG -s 2 -pf y ) \
            > "$dir/antechamber.log" 2>&1 &
        progress $!
        wait $! || true

        # Some antechamber builds do not accept sdf; fall back to mol2.
        if [[ ! -s "$dir/lig.mol2" ]]; then
            echo "sdf input rejected, retrying with mol2..."
            ( cd "$dir" && antechamber -i lig_in.mol2 -fi mol2 -o lig.mol2 -fo mol2 \
                -c bcc -nc "$nc" -at gaff2 -rn LIG -s 2 -pf y ) \
                >> "$dir/antechamber.log" 2>&1 || true
        fi
        done_
        local t_charge=$LAST_STEP_SECONDS

        [[ -s "$dir/lig.mol2" ]] || {
            echo "ERROR: antechamber failed for $cpd. See $dir/antechamber.log"
            tail -n 15 "$dir/antechamber.log"
            return 1
        }
        ( cd "$dir" && parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2 ) \
            > "$dir/parmchk2.log" 2>&1
    fi
    local t_charge="${t_charge:-0}"

    # ---- topologies: complex / receptor / ligand ----
    cat > "$dir/tleap.in" <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff2

loadamberparams ../FMN.frcmod
loadamberparams lig.frcmod

FMN = loadmol2 ../FMN.mol2
LIG = loadmol2 lig.mol2

prot = loadpdb ../protein_amber.pdb

rec = combine { prot FMN }
com = combine { prot FMN LIG }

set default PBRadii mbondi2

saveamberparm LIG ligand.prmtop ligand.inpcrd
saveamberparm rec receptor.prmtop receptor.inpcrd
saveamberparm com complex.prmtop complex.inpcrd
savepdb com complex_start.pdb

quit
EOF

    begin "tleap topologies"
    ( cd "$dir" && tleap -f tleap.in ) > "$dir/tleap.log" 2>&1 || true
    done_
    local t_tleap=$LAST_STEP_SECONDS

    for f in complex.prmtop receptor.prmtop ligand.prmtop complex.inpcrd; do
        [[ -s "$dir/$f" ]] || {
            echo "ERROR: tleap did not create $f for $cpd."
            grep -Ei 'error|fatal|could not|unknown' "$dir/tleap.log" | head -n 20
            return 1
        }
    done

    # ---- minimisation (implicit solvent) ----
    cat > "$dir/min.in" <<EOF
Minimisation of the docked complex, GB implicit solvent
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
    local t_min=$LAST_STEP_SECONDS

    [[ -s "$dir/complex_min.rst7" ]] || {
        echo "ERROR: minimisation failed for $cpd."
        tail -n 20 "$dir/min.out" 2>/dev/null || cat "$dir/sander.log"
        return 1
    }

    # ---- single-frame trajectory for MMPBSA.py ----
    cat > "$dir/traj.in" <<'EOF'
trajin complex_min.rst7
trajout complex_min.nc netcdf
go
quit
EOF
    ( cd "$dir" && cpptraj -p complex.prmtop -i traj.in ) \
        > "$dir/cpptraj.log" 2>&1

    [[ -s "$dir/complex_min.nc" ]] || {
        echo "ERROR: cpptraj could not write the trajectory for $cpd."
        tail -n 10 "$dir/cpptraj.log"
        return 1
    }

    # ---- MM/GBSA ----
    cat > "$dir/mmgbsa.in" <<EOF
MM/GBSA on the minimised docked pose
&general
  startframe = 1, endframe = 1, interval = 1,
  verbose = 2, keep_files = 0,
/
&gb
  igb = ${IGB}, saltcon = ${SALT},
/
EOF

    begin "MM/GBSA"
    ( cd "$dir" && MMPBSA.py -O -i mmgbsa.in -o mmgbsa.dat -eo mmgbsa_frames.csv \
        -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
        -y complex_min.nc ) > "$dir/mmpbsa.log" 2>&1 &
    progress $!
    wait $! || true
    done_
    local t_gbsa=$LAST_STEP_SECONDS

    [[ -s "$dir/mmgbsa.dat" ]] || {
        echo "ERROR: MMPBSA.py failed for $cpd."
        tail -n 25 "$dir/mmpbsa.log"
        return 1
    }

    ( cd "$dir" && cpptraj -p complex.prmtop <<'EOF'
trajin complex_min.rst7
trajout complex_min.pdb pdb
go
quit
EOF
    ) > "$dir/cpptraj_pdb.log" 2>&1 || true

    stamp_pdb "$dir/complex_min.pdb" "${cpd} minimised complex (MM/GBSA, ${RUN} pose)"
    stamp_pdb "$dir/complex_start.pdb" "${cpd} docked complex before minimisation"

    local t_total=$(( SECONDS - t_start ))
    printf '  %-34s%s\n' "TOTAL for ${cpd}" "$(hms "$t_total")"

    {
        echo "compound ${cpd}"
        echo "  AM1-BCC charges  $(hms "${t_charge:-0}")"
        echo "  tleap            $(hms "${t_tleap:-0}")"
        echo "  minimisation     $(hms "${t_min:-0}") (${MAXCYC} cycles)"
        echo "  MM/GBSA          $(hms "${t_gbsa:-0}")"
        echo "  total            $(hms "$t_total")"
    } >> "$WORK/timings.txt"
}

RUN_START=$SECONDS
: > "$WORK/timings.txt"

FAILED=()
for cpd in "${COMPOUNDS[@]}"; do
    run_compound "$cpd" || FAILED+=("$cpd")
done

TOTAL_SECONDS=$(( SECONDS - RUN_START ))

# ============================================================
# 4D. Collect, compare, graph
# ============================================================

python3 - "$WORK" "$TXT_DIR" "$POSE_PREFIX" "$RUN" "${COMPOUNDS[@]}" <<'PY'
import math
import re
import sys
from datetime import datetime
from pathlib import Path

work, txt_dir, prefix, run = sys.argv[1:5]
cpds = sys.argv[5:]
work, txt_dir = Path(work), Path(txt_dir)

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
NOW = f"{datetime.now():%Y-%m-%d %H:%M:%S}"
HEAD = ["=" * 72, *IDENTITY, f"Date: {NOW}", "=" * 72, ""]
FOOTER = ("siam chowdhury | Computational and Medicinal Chemistry | "
          "Dr. Alam's Research Team (www.alamresearch.org) | Arkansas State University")

TERMS = [
    ("VDWAALS", "van der Waals"),
    ("EEL", "electrostatic"),
    ("EGB", "polar solvation"),
    ("ESURF", "nonpolar solvation"),
    ("DELTA G gas", "gas-phase total"),
    ("DELTA G solv", "solvation total"),
    ("DELTA TOTAL", "MM/GBSA binding energy"),
]


def parse_mmgbsa(path: Path):
    """Read the DELTAS section of an MMPBSA.py output file."""
    text = path.read_text(errors="replace")
    block = text.split("Differences (Complex - Receptor - Ligand)")[-1]
    out = {}
    for key, _ in TERMS:
        m = re.search(
            rf"^{re.escape(key)}\s+(-?\d+\.\d+)", block, re.M
        )
        if m:
            out[key] = float(m.group(1))
    return out


rows = []
for c in cpds:
    dat = work / c / "mmgbsa.dat"
    if not dat.is_file():
        continue
    vals = parse_mmgbsa(dat)
    if "DELTA TOTAL" not in vals:
        continue

    vina, heavy, contacts, shared, site, fmn, polar = (None,) * 7
    tsv = Path(f"{prefix}{c}/{c}_summary.tsv")
    if tsv.is_file():
        v = tsv.read_text().split("\t")
        vina, heavy = float(v[1]), int(v[2])
        contacts, shared, site = int(v[6]), int(v[7]), int(v[8])
        fmn, polar = v[9], int(v[10].strip())

    rows.append({
        "cpd": c, "vals": vals, "vina": vina, "heavy": heavy,
        "contacts": contacts, "shared": shared, "site": site,
        "fmn": fmn, "polar": polar,
    })

if not rows:
    raise SystemExit("ERROR: no MM/GBSA results could be read.")

rows.sort(key=lambda r: r["vals"]["DELTA TOTAL"])

lines = HEAD + [
    f"36 MM/GBSA RESULTS: 7L00 FabK ({run} poses)",
    "-" * 72,
    "Single minimised structure, GB implicit solvent (igb=5, 0.15 M salt).",
    "Entropy is NOT included, so these are relative, not absolute, energies.",
    "",
    f"{'RANK':<5}{'COMPOUND':<11}{'MMGBSA':>9}{'VDW':>9}{'EEL':>10}"
    f"{'EGB':>10}{'ESURF':>9}{'VINA':>8}{'LE_GBSA':>9}",
]
for i, r in enumerate(rows, start=1):
    v = r["vals"]
    le = -v["DELTA TOTAL"] / r["heavy"] if r["heavy"] else float("nan")
    lines.append(
        f"{i:<5}{r['cpd']:<11}{v['DELTA TOTAL']:>9.2f}{v.get('VDWAALS', 0):>9.2f}"
        f"{v.get('EEL', 0):>10.2f}{v.get('EGB', 0):>10.2f}{v.get('ESURF', 0):>9.2f}"
        f"{(r['vina'] if r['vina'] is not None else float('nan')):>8.2f}{le:>9.3f}"
    )

lines += [
    "",
    "TERMS (kcal/mol, all are complex minus receptor minus ligand)",
    "  VDW    van der Waals: shape fit and buried surface",
    "  EEL    electrostatic: charge-charge, salt bridges, H-bonds",
    "  EGB    polar solvation: cost of desolvating charges on binding",
    "         (opposes EEL; a strong salt bridge pays a large EGB penalty)",
    "  ESURF  nonpolar solvation: burial of hydrophobic surface",
    "  MMGBSA total binding energy, more negative = stronger predicted binding",
    "  LE_GBSA  -MMGBSA / heavy atoms",
    "",
    "CROSS-METHOD COMPARISON",
    "-" * 72,
    f"{'COMPOUND':<11}{'MMGBSA':>9}{'RANK':>6}{'VINA':>8}{'RANK':>6}"
    f"{'FMN':>6}{'POLAR':>7}{'SITE':>8}",
]
vina_order = sorted(
    [r for r in rows if r["vina"] is not None], key=lambda r: r["vina"]
)
vina_rank = {r["cpd"]: i for i, r in enumerate(vina_order, start=1)}
for i, r in enumerate(rows, start=1):
    lines.append(
        f"{r['cpd']:<11}{r['vals']['DELTA TOTAL']:>9.2f}{i:>6}"
        f"{(r['vina'] if r['vina'] is not None else float('nan')):>8.2f}"
        f"{vina_rank.get(r['cpd'], 0):>6}{str(r['fmn']):>6}{r['polar']:>7}"
        f"{str(r['shared']) + '/' + str(r['site']):>8}"
    )

agree = [r["cpd"] for i, r in enumerate(rows, start=1) if vina_rank.get(r["cpd"]) == i]
lines += [
    "",
    f"Methods agree on the rank of: {', '.join(agree) if agree else 'no compound'}",
    "",
    "Where MM/GBSA and Vina disagree, MM/GBSA is usually the better guide,",
    "because it treats electrostatics and desolvation explicitly, while the",
    "Vina score weights buried hydrophobic surface heavily. Neither includes",
    "entropy, and neither is an experimental binding free energy.",
    "",
]

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

pairs = [(r["vals"]["DELTA TOTAL"], mic[r["cpd"]]) for r in rows if r["cpd"] in mic]
if len(pairs) >= 3:
    ys = [p[0] for p in pairs]
    xs = [math.log10(p[1]) for p in pairs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    rv = sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")
    lines += [
        f"MIC vs MM/GBSA (n={len(pairs)}): Pearson r = {rv:.3f} against log10(MIC)",
        "Expected: lower MIC with more negative MM/GBSA, so r > 0.",
        "",
    ]

(txt_dir / "36_mmgbsa_results.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines[len(HEAD):]))

# ---- per-compound detail files ----
for r in rows:
    v = r["vals"]
    detail = HEAD + [
        f"37 MM/GBSA DETAIL: {r['cpd']} ({run} pose)",
        "-" * 72,
        f"{'TERM':<22}{'kcal/mol':>12}   MEANING",
    ]
    for key, meaning in TERMS:
        if key in v:
            detail.append(f"{key:<22}{v[key]:>12.2f}   {meaning}")
    detail += [
        "",
        f"Heavy atoms            {r['heavy']}",
        f"Vina score             {r['vina']:.2f} kcal/mol" if r["vina"] is not None else "",
        f"FMN contact (pose 1)   {r['fmn']}",
        f"Polar contacts         {r['polar']}",
        f"XCJ pocket residues    {r['shared']} of {r['site']}",
        "",
        "FILES",
        f"  {work}/{r['cpd']}/complex_min.pdb    minimised complex",
        f"  {work}/{r['cpd']}/mmgbsa.dat         full MMPBSA.py output",
        "",
        "No entropy term is included, so this is a relative ranking energy,",
        "not an experimental binding free energy.",
    ]
    (txt_dir / f"37_{r['cpd']}_mmgbsa_detail.txt").write_text(
        "\n".join(l for l in detail if l != "") + "\n"
    )

# ---- graphs ----
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit(0)

graph_dir = Path("graphs" if run == "global" else f"graphs_{run}")
graph_dir.mkdir(exist_ok=True)


def finish(fig, path):
    fig.text(0.5, 0.005, FOOTER, ha="center", va="bottom", fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


names = [r["cpd"] for r in rows]
totals = [r["vals"]["DELTA TOTAL"] for r in rows]

fig = plt.figure(figsize=(max(6, 1.4 * len(names) + 3), 6))
bars = plt.bar(names, totals, color="#55A868")
for b, s in zip(bars, totals):
    plt.annotate(f"{s:.1f}", (b.get_x() + b.get_width() / 2, s), xytext=(0, 4),
                 textcoords="offset points", ha="center", va="bottom", fontsize=9)
plt.ylabel("MM/GBSA binding energy (kcal/mol)")
plt.title(f"MM/GBSA Binding Energy: 7L00 FabK [{run} poses]")
plt.ylim(0, min(totals) * 1.15)
plt.grid(True, axis="y", linestyle="--", alpha=0.4)
finish(fig, graph_dir / "mmgbsa_binding_energy.png")

keys = [k for k, _ in TERMS if k in {"VDWAALS", "EEL", "EGB", "ESURF"}]
width = 0.8 / len(keys)
fig = plt.figure(figsize=(max(7, 1.8 * len(names) + 4), 6))
for j, k in enumerate(keys):
    plt.bar(
        [i + j * width - 0.4 for i in range(len(rows))],
        [r["vals"].get(k, 0) for r in rows],
        width=width, label=k,
    )
plt.xticks(range(len(rows)), names)
plt.ylabel("Energy component (kcal/mol)")
plt.title(f"MM/GBSA Energy Components: 7L00 FabK [{run} poses]")
plt.axhline(0, color="black", linewidth=0.8)
plt.legend()
plt.grid(True, axis="y", linestyle="--", alpha=0.4)
finish(fig, graph_dir / "mmgbsa_components.png")

if len(rows) >= 2 and all(r["vina"] is not None for r in rows):
    fig = plt.figure(figsize=(7, 6))
    for r in rows:
        plt.scatter(r["vina"], r["vals"]["DELTA TOTAL"], s=70)
        plt.annotate(r["cpd"], (r["vina"], r["vals"]["DELTA TOTAL"]),
                     xytext=(6, 4), textcoords="offset points", fontsize=9)
    plt.xlabel("Vina score (kcal/mol)")
    plt.ylabel("MM/GBSA binding energy (kcal/mol)")
    plt.title(f"Vina vs MM/GBSA: 7L00 FabK [{run} poses]")
    plt.grid(True, linestyle="--", alpha=0.4)
    finish(fig, graph_dir / "vina_vs_mmgbsa.png")

if len(pairs) >= 2:
    fig = plt.figure(figsize=(7, 6))
    for r in rows:
        if r["cpd"] in mic:
            plt.scatter(mic[r["cpd"]], r["vals"]["DELTA TOTAL"], s=70)
            plt.annotate(r["cpd"], (mic[r["cpd"]], r["vals"]["DELTA TOTAL"]),
                         xytext=(6, 4), textcoords="offset points", fontsize=9)
    plt.xscale("log")
    plt.xlabel("MIC (ug/mL, log scale)")
    plt.ylabel("MM/GBSA binding energy (kcal/mol)")
    plt.title("MIC vs MM/GBSA: 7L00 FabK")
    plt.grid(True, linestyle="--", alpha=0.4)
    finish(fig, graph_dir / "MIC_vs_mmgbsa.png")

print(f"\nGraphs written to {graph_dir}/")
PY

echo
echo "============================================================"
echo "STEP 4 COMPLETED"
echo "============================================================"
echo "Total run time: $(hms "${TOTAL_SECONDS:-0}")"
echo
if [[ -s "$WORK/timings.txt" ]]; then
    {
        header_txt
        echo "38 RUN TIMINGS"
        echo "------------------------------------------------------------"
        echo "Settings: MAXCYC=${MAXCYC}, NCYC=${NCYC}, igb=${IGB}, chains ${CHAINS}"
        echo
        cat "$WORK/timings.txt"
        echo
        echo "total run     $(hms "${TOTAL_SECONDS:-0}")"
    } > "$TXT_DIR/38_run_timings.txt"
    cat "$WORK/timings.txt"
    echo
fi
if (( ${#FAILED[@]} )); then
    echo "FAILED: ${FAILED[*]}   (see ${WORK}/<compound>/*.log)"
fi
echo "Results : cat $TXT_DIR/36_mmgbsa_results.txt"
echo "Detail  : cat $TXT_DIR/37_<compound>_mmgbsa_detail.txt"
echo "Setup   : cat $TXT_DIR/35_mmgbsa_receptor_setup.txt"
echo "Structures: ${WORK}/<compound>/complex_min.pdb"
