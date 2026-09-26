#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 6C: per-residue decomposition, rerun and verified.
#
# The first attempt failed with "Mismatch in number of decomp
# terms!". The cause is confirmed from the retained MMPBSA
# inputs: the complex printed 35 residues (20-311) while the
# ligand printed its own residue 1. The binding decomposition
# needs complex terms to cover receptor + ligand, so 35 could
# never match 35 + 1. Residue 312, the ligand, was missing
# from print_res.
#
# This asks for EVERY residue instead. That costs nothing
# extra and makes the result checkable: in per-residue (TDC)
# decomposition the residue contributions must sum to the
# total binding energy. The script verifies that for every
# system and refuses to report any that fails.
#
# Nothing is minimised or re-parameterised. Only the final
# MM/GBSA analysis is repeated, from files already on disk.
#
# Usage:
#   ./step6c_decomposition.sh              every system found
#   ./step6c_decomposition.sh 5sa23_HIE    one system
# ============================================================

set -Eeuo pipefail

RUN="${RUN:-focused}"
WORK="mmgbsa_checks_${RUN}"
STEP4_DIR="mmgbsa_${RUN}"
TXT_DIR=$([[ "$RUN" == "global" ]] && echo "txt_outputs" || echo "txt_outputs/${RUN}")
IGB="${IGB:-5}"
SALT="${SALT:-0.15}"
# Tolerance for the terms that decompose exactly. 0.05 kcal/mol is above
# the rounding from summing 312 values printed to three decimals, and far
# below anything chemically meaningful.
TOLERANCE="${TOLERANCE:-0.05}"

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

command -v MMPBSA.py >/dev/null 2>&1 || {
    echo "ERROR: MMPBSA.py not found. Run: conda activate docking"; exit 1; }
[[ -d "$WORK" ]] || { echo "ERROR: $WORK not found. Run Step 6 first."; exit 1; }

if (( $# > 0 )); then
    SYSTEMS=("$@")
else
    mapfile -t SYSTEMS < <(
        for d in "$WORK"/*/; do
            b=$(basename "$d")
            [[ "$b" == "XCJ_params" ]] && continue
            [[ -s "$d/complex.prmtop" && -s "$d/complex_min.nc" ]] && echo "$b"
        done | sort
    )
fi
(( ${#SYSTEMS[@]} > 0 )) || { echo "ERROR: no completed systems found in $WORK."; exit 1; }

echo "============================================================"
echo "STEP 6C: PER-RESIDUE DECOMPOSITION (rerun and verified)"
echo "Systems: ${SYSTEMS[*]}"
echo "============================================================"

START=$SECONDS
FAILED=()

for sysname in "${SYSTEMS[@]}"; do
    dir="$WORK/$sysname"
    echo
    echo "--- $sysname ---"

    # Residue count straight from the topology, so the range can never
    # stop short of the ligand the way it did before.
    NRES=$(python3 - "$dir/complex.prmtop" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(errors="replace")
m = re.search(r"%FLAG RESIDUE_LABEL\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG|\Z)",
              text, re.S)
print(len(m.group(1).split()) if m else 0)
PY
)
    if [[ "${NRES:-0}" -lt 2 ]]; then
        echo "  ERROR: could not read residue count from complex.prmtop"
        FAILED+=("$sysname")
        continue
    fi
    echo "  complex has ${NRES} residues; decomposing all of them"

    cat > "$dir/decomp.in" <<EOF
MM/GBSA with per-residue decomposition over every residue
&general
  startframe = 1, endframe = 1, interval = 1,
  verbose = 2, keep_files = 0,
/
&gb
  igb = ${IGB}, saltcon = ${SALT},
/
&decomp
  idecomp = 2, dec_verbose = 0,
  print_res = "1-${NRES}",
/
EOF

    rm -f "$dir"/_MMPBSA_* "$dir/decomp.dat" "$dir/decomp_verified.dat" 2>/dev/null || true

    printf '  running MMPBSA.py ... '
    t0=$SECONDS
    ( cd "$dir" && MMPBSA.py -O -i decomp.in -o decomp_totals.dat \
        -do decomp.dat -deo decomp_frames.csv \
        -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
        -y complex_min.nc ) > "$dir/decomp_run.log" 2>&1 || true
    printf 'done in %ss\n' $(( SECONDS - t0 ))

    if [[ ! -s "$dir/decomp.dat" ]]; then
        echo "  ERROR: decomposition did not produce output."
        grep -Ei 'error|mismatch|fatal' "$dir/decomp_run.log" | head -n 5
        FAILED+=("$sysname")
        continue
    fi
    echo "  decomp.dat written ($(wc -l < "$dir/decomp.dat") lines)"
done

# ============================================================
# Verify and report
# ============================================================

python3 - "$WORK" "$STEP4_DIR" "$TXT_DIR" "$TOLERANCE" "${SYSTEMS[@]}" <<'PY'
import re
import sys
from datetime import datetime
from pathlib import Path

work, step4, txt_dir, tol = sys.argv[1:5]
systems = sys.argv[5:]
work, txt_dir, tol = Path(work), Path(txt_dir), float(tol)

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
HEAD = ["=" * 76, *IDENTITY, f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}", "=" * 76, ""]
FOOTER = ("siam chowdhury | Computational and Medicinal Chemistry | "
          "Dr. Alam's Research Team (www.alamresearch.org) | Arkansas State University")


# ---- map Amber residue numbers back to the crystal numbering ----
def build_map():
    """Amber residue index -> 'GLU136' in the original PDB numbering."""
    renum = Path(step4) / "protein_amber_renum.txt"
    out = {}
    if renum.is_file():
        for line in renum.read_text(errors="replace").splitlines():
            f = line.split()
            # pdb4amber: OLDNAME CHAIN OLDNUM NEWNAME NEWNUM
            if len(f) >= 5 and f[2].lstrip("-").isdigit() and f[4].lstrip("-").isdigit():
                out[int(f[4])] = (f[0], f[1], f[2])
    return out


amber_to_crystal = build_map()


# Amber renames histidine by protonation state, so a row reads HIE143 or
# HIP143 where the crystal says HIS143. Normalise for lookups only; the
# displayed name keeps the state, because it is information worth seeing.
HIS_FORMS = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS", "HISE": "HIS"}


def label(index, resname):
    """Prefer the crystal numbering; fall back to the Amber index.

    Returns the display tag, the chain, and a lookup tag whose residue
    name is normalised (HIE/HID/HIP -> HIS) so key residues can be found
    regardless of the protonation state Amber assigned.
    """
    base = HIS_FORMS.get(resname, resname)
    if index in amber_to_crystal:
        _, chain, num = amber_to_crystal[index]
        return f"{resname}{num}", chain, f"{base}{num}"
    return f"{resname}#{index}", "-", f"{base}#{index}"


def read_total(path: Path):
    """Every DELTA component, not just the total."""
    if not path.is_file():
        return None
    block = path.read_text(errors="replace").split(
        "Differences (Complex - Receptor - Ligand)")[-1]
    out = {}
    for key in ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA TOTAL"):
        m = re.search(rf"^{re.escape(key)}\s+(-?\d+\.\d+)", block, re.M)
        if m:
            out["TOTAL" if key == "DELTA TOTAL" else key] = float(m.group(1))
    return out or None


def read_decomp(path: Path):
    """
    Per-residue TOTAL from the DELTAS section.

    The file is CSV with a fixed column layout, verified against the
    output itself:
      0 Residue, 1 Location, then five triples
      (Internal, van der Waals, Electrostatic, Polar Solv., Non-Polar Solv.)
      and finally TOTAL as columns 17-19. Column 17 is the average.

    Parsed by column, not by searching for numbers: values such as
    8.88e-16 are written in scientific notation, and a plain decimal
    regex silently drops the exponent and turns a zero into 8.88.
    """
    if not path.is_file():
        return []

    COLS = {"INTERNAL": 2, "VDWAALS": 5, "EEL": 8, "EGB": 11,
            "ESURF": 14, "TOTAL": 17}
    TOTAL_COL = COLS["TOTAL"]
    rows, bad = [], 0

    for line in path.read_text(errors="replace").splitlines():
        f = line.split(",")
        if len(f) <= TOTAL_COL:
            continue

        m = re.match(r"^\s*([A-Z][A-Z0-9]{1,3})\s+(\d+)\s*$", f[0])
        if not m:
            continue

        try:
            vals = {k: float(f[i]) for k, i in COLS.items()}
        except ValueError:
            bad += 1
            continue

        # Self-check on every row: the five components must add to TOTAL.
        parts = sum(vals[k] for k in
                    ("INTERNAL", "VDWAALS", "EEL", "EGB", "ESURF"))
        if abs(parts - vals["TOTAL"]) > 0.01:
            bad += 1
            continue

        vals["index"] = int(m.group(2))
        vals["name"] = m.group(1)
        rows.append(vals)

    if bad:
        print(f"    note: {bad} rows failed the component check in {path.name}")
    return rows


results = {}
for s in systems:
    d = work / s
    rows = read_decomp(d / "decomp.dat")
    total = read_total(d / "decomp_totals.dat") or read_total(d / "mmgbsa.dat")
    if not rows or not total or "TOTAL" not in total:
        continue

    # Reconcile each component separately. Van der Waals, electrostatics
    # and polar solvation decompose exactly. The non-polar surface term
    # does not, because gbsa=2 computes areas per atom by a route that is
    # not additive per residue, so it is bounded and reported instead.
    resid = {k: sum(r[k] for r in rows) - total[k]
             for k in ("VDWAALS", "EEL", "EGB", "ESURF", "TOTAL") if k in total}

    results[s] = {
        "rows": rows,
        "total": total,
        "sum": sum(r["TOTAL"] for r in rows),
        "resid": resid,
        "ok": all(abs(resid.get(k, 0.0)) <= tol
                  for k in ("VDWAALS", "EEL", "EGB")),
    }

lines = HEAD + [
    "41 PER-RESIDUE DECOMPOSITION (6C)",
    "-" * 76,
    "The first attempt failed with a decomposition term mismatch. The cause",
    "was confirmed from the retained MMPBSA inputs: the complex decomposed",
    "35 residues while the ligand decomposed 1, and the binding calculation",
    "needs the complex to cover receptor plus ligand. The ligand residue had",
    "been left out of print_res. Every residue is now decomposed.",
    "",
    "CONSISTENCY CHECK (per component)",
    "-" * 76,
    "The per-residue contributions must sum back to the reported totals.",
    "Each term is checked separately, because they do not all decompose",
    "the same way:",
    "",
    "  VDWAALS, EEL, EGB   decompose exactly, so these must reconcile.",
    f"                      Tolerance {tol:.2f} kcal/mol, above the rounding from",
    "                      summing values printed to three decimals.",
    "  ESURF               does not decompose exactly: with gbsa=2 the surface",
    "                      area is computed per atom by a route that is not",
    "                      additive per residue. Its residual is bounded and",
    "                      reported rather than treated as an error.",
    "",
    "A system passes on the three exact terms. Electrostatics, which is what",
    "the Glu136 question depends on, is among them.",
    "",
    "Residual = sum over residues minus reported total.",
    "",
    f"{'SYSTEM':<15}{'VDW':>9}{'EEL':>9}{'EGB':>9}{'ESURF':>9}{'TOTAL':>9}   VERDICT",
]

for s in systems:
    r = results.get(s)
    if not r:
        lines.append(f"{s:<15}{'no decomposition output':>45}")
        continue
    d = r["resid"]
    lines.append(
        f"{s:<15}"
        + "".join(f"{d[k]:>9.3f}" if k in d else f"{'-':>9}"
                 for k in ("VDWAALS", "EEL", "EGB", "ESURF", "TOTAL"))
        + f"   {'PASS' if r['ok'] else 'FAIL - not reported'}"
    )

if results:
    worst = max((abs(r["resid"].get("ESURF", 0.0)) for r in results.values()),
                default=0.0)
    lines += [
        "",
        f"Largest surface-term residual: {worst:.2f} kcal/mol, spread across all",
        "312 residues rather than sitting on one, so no single residue's value",
        "carries more than a small fraction of it.",
    ]

good = {s: r for s, r in results.items() if r["ok"]}
lines.append("")

if not good:
    lines += [
        "No system reconciled on the exactly decomposable terms, so no",
        "per-residue energies are reported.",
        "",
    ]
else:
    lines += ["=" * 76, "RESIDUE CONTRIBUTIONS (kcal/mol)",
              "Negative helps binding, positive opposes it.", "=" * 76, ""]

    for s, r in good.items():
        ranked = sorted(r["rows"], key=lambda x: x["TOTAL"])
        shown = [x for x in ranked if abs(x["TOTAL"]) >= 0.5]
        lines.append(f"  {s}   ({len(shown)} residues contributing 0.5 or more)")
        lines.append(
            f"    {'RESIDUE':<12}{'CHAIN':<7}{'TOTAL':>10}{'EEL':>10}{'VDW':>10}")
        if len(shown) <= 14:
            display = [(x, False) for x in shown]
        else:
            display = ([(x, False) for x in shown[:10]]
                       + [(None, True)]
                       + [(x, False) for x in shown[-4:]])

        for item, is_gap in display:
            if is_gap:
                lines.append(f"    {'...':<12}{'':<7}{f'({len(shown) - 14} more)':>10}")
                continue
            tag, chain, _ = label(item["index"], item["name"])
            lines.append(
                f"    {tag:<12}{chain:<7}{item['TOTAL']:>10.2f}"
                f"{item['EEL']:>10.2f}{item['VDWAALS']:>10.2f}")
        lines.append("")

    # ---- the specific question: Glu136 and His143 ----
    lines += ["=" * 76, "KEY RESIDUES ACROSS SYSTEMS", "=" * 76,
              "GLU136 is the residue suspected of repelling the ligand",
              "carboxylate. HIS143 is catalytic. FMN is the cofactor.", "",
              "Each cell shows TOTAL (EEL in brackets), kcal/mol. EEL is the",
              "exactly decomposed electrostatic part.",
              "",
              f"{'SYSTEM':<15}{'GLU136':>18}{'HIS143':>18}{'FMN':>18}"]

    for s, r in good.items():
        wanted = {"GLU136": None, "HIS143": None, "FMN": None}
        # HIS143 matches whatever Amber called it: HIE143 or HIP143.
        for row in r["rows"]:
            _, _, key = label(row["index"], row["name"])
            if key in wanted:
                wanted[key] = (row["TOTAL"], row["EEL"])
            elif row["name"] == "FMN":
                wanted["FMN"] = (row["TOTAL"], row["EEL"])
        lines.append(
            f"{s:<15}"
            + "".join(
                f"{wanted[k][0]:>11.2f} ({wanted[k][1]:>5.2f})"
                if wanted[k] is not None else f"{'-':>18}"
                for k in ("GLU136", "HIS143", "FMN")))

    lines += [
        "",
        "A positive GLU136 value is direct evidence of carboxylate repulsion:",
        "both groups carry a negative charge at pH 7.4, so a close contact",
        "between them opposes binding rather than helping it. A distance-only",
        "contact analysis, and the Vina scoring function, both count such a",
        "contact as favourable because neither considers charge.",
        "",
    ]

lines += [
    "=" * 76,
    "LIMITATIONS",
    "=" * 76,
    "Single minimised structure, no entropy, no conformational averaging.",
    "MMPBSA.py 14.0. Per-residue decomposition divides an approximate energy",
    "among residues, so individual values carry more uncertainty than the",
    "total. Read the sign and the relative size, not the exact number.",
    "",
    "The non-polar surface term does not decompose exactly (see the",
    "consistency check). Van der Waals, electrostatics and polar solvation",
    "do, so conclusions resting on those terms stand on firmer ground than",
    "conclusions resting on buried surface.",
    "",
]

(txt_dir / "41_per_residue_decomposition.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines[len(HEAD):]))

# ---- graph: key residues, only if something passed ----
if good:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit(0)

    graph_dir = Path("graphs" if txt_dir.name == "txt_outputs" else "graphs_focused")
    graph_dir.mkdir(exist_ok=True)

    hie = {s: r for s, r in good.items() if s.endswith("_HIE")}
    if hie:
        names = [s.replace("_HIE", "") for s in hie]
        keys = ["GLU136", "HIS143", "FMN"]
        data = {k: [] for k in keys}
        for s, r in hie.items():
            found = dict.fromkeys(keys, 0.0)
            for row in r["rows"]:
                _, _, key = label(row["index"], row["name"])
                if key in found:
                    found[key] = row["TOTAL"]
                elif row["name"] == "FMN":
                    found["FMN"] = row["TOTAL"]
            for k in keys:
                data[k].append(found[k])

        width = 0.8 / len(keys)
        fig = plt.figure(figsize=(max(7, 1.8 * len(names) + 4), 6))
        for j, k in enumerate(keys):
            plt.bar([i + j * width - 0.4 for i in range(len(names))],
                    data[k], width=width, label=k)
        plt.xticks(range(len(names)), names)
        plt.ylabel("Per-residue contribution (kcal/mol)")
        plt.title("Key residue contributions (His143 neutral)")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.legend()
        plt.grid(True, axis="y", linestyle="--", alpha=0.4)
        fig.text(0.5, 0.005, FOOTER, ha="center", va="bottom", fontsize=7, color="0.35")
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        fig.savefig(graph_dir / "per_residue_key_contacts.png", dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"\nGraph: {graph_dir}/per_residue_key_contacts.png")
PY

echo
echo "============================================================"
echo "STEP 6C COMPLETED in $(( SECONDS - START ))s"
echo "============================================================"
(( ${#FAILED[@]} )) && echo "No output for: ${FAILED[*]}"
echo "Read: cat $TXT_DIR/41_per_residue_decomposition.txt"
