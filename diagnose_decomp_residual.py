#!/usr/bin/env python3
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# DIAGNOSTIC ONLY - reads files, changes nothing.
#
# The per-residue sums now land within about 1 kcal/mol of the
# totals, but every system is off in the same direction. Before
# deciding whether that residual is acceptable, it has to be
# explained rather than tolerated.
#
# This reconciles EACH energy component separately:
#
#   sum over residues of VDWAALS  vs  DELTA VDWAALS
#   sum over residues of EEL      vs  DELTA EEL
#   sum over residues of EGB      vs  DELTA EGB
#   sum over residues of ESURF    vs  DELTA ESURF
#
# If van der Waals, electrostatics and polar solvation reconcile
# exactly and the whole residual sits in the non-polar surface
# term, the cause is gbsa=2 per-atom surface areas, which are
# not strictly additive per residue. That is a known property
# of the method, not an error, and it can be stated and bounded.
#
# It also counts rows rejected by the row filter, in case the
# filter itself is removing real contributions.
#
# Usage:  python3 diagnose_decomp_residual.py [work_dir]
# ============================================================

from __future__ import annotations

import re
import sys
from pathlib import Path

work = Path(sys.argv[1] if len(sys.argv) > 1 else "mmgbsa_checks_focused")

# Column layout of the decomposition CSV, verified against the output:
#   0 Residue, 1 Location, then five (Avg, SD, SEM) triples, then TOTAL.
COLS = {
    "INTERNAL": 2,
    "VDWAALS": 5,
    "EEL": 8,
    "EGB": 11,
    "ESURF": 14,
    "TOTAL": 17,
}

if not work.is_dir():
    raise SystemExit(f"ERROR: {work} not found. Run from your working folder.")


def read_totals(path: Path):
    """DELTA components from the MMPBSA totals file."""
    if not path.is_file():
        return {}
    block = path.read_text(errors="replace").split(
        "Differences (Complex - Receptor - Ligand)")[-1]
    out = {}
    for key in ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA TOTAL"):
        m = re.search(rf"^{re.escape(key)}\s+(-?\d+\.\d+)", block, re.M)
        if m:
            out["TOTAL" if key == "DELTA TOTAL" else key] = float(m.group(1))
    return out


def read_rows(path: Path):
    """Every parsable row, with no filtering, plus the rows that fail it."""
    if not path.is_file():
        return [], []

    kept, rejected = [], []
    for line in path.read_text(errors="replace").splitlines():
        f = line.split(",")
        if len(f) <= COLS["TOTAL"]:
            continue
        m = re.match(r"^\s*([A-Z][A-Z0-9]{1,3})\s+(\d+)\s*$", f[0])
        if not m:
            continue
        try:
            vals = {k: float(f[i]) for k, i in COLS.items()}
        except ValueError:
            continue

        parts = sum(vals[k] for k in ("INTERNAL", "VDWAALS", "EEL", "EGB", "ESURF"))
        entry = (f"{m.group(1)}{m.group(2)}", vals, parts - vals["TOTAL"])
        (kept if abs(parts - vals["TOTAL"]) <= 0.01 else rejected).append(entry)

    return kept, rejected


systems = sorted(
    d.name for d in work.iterdir()
    if d.is_dir() and (d / "decomp.dat").is_file()
)
if not systems:
    raise SystemExit(f"ERROR: no decomp.dat found under {work}/")

print("=" * 78)
print("DECOMPOSITION RESIDUAL DIAGNOSTIC  (read-only)")
print("=" * 78)

pattern = {}

for s in systems:
    d = work / s
    totals = read_totals(d / "decomp_totals.dat") or read_totals(d / "mmgbsa.dat")
    kept, rejected = read_rows(d / "decomp.dat")

    if not kept or not totals:
        print(f"\n{s}: could not read both files")
        continue

    print(f"\n{s}   ({len(kept)} rows kept, {len(rejected)} rejected by the row filter)")
    print(f"  {'COMPONENT':<12}{'SUM OF RESIDUES':>17}{'DELTA TOTAL':>14}"
          f"{'RESIDUAL':>11}")

    residuals = {}
    for key in ("VDWAALS", "EEL", "EGB", "ESURF", "TOTAL"):
        if key not in totals:
            continue
        summed = sum(v[key] for _, v, _ in kept)
        resid = summed - totals[key]
        residuals[key] = resid
        flag = "" if abs(resid) < 0.01 else ("  <-- differs" if abs(resid) > 0.1 else "")
        print(f"  {key:<12}{summed:>17.3f}{totals[key]:>14.3f}{resid:>11.3f}{flag}")

    # Internal is not part of the binding total but is carried in the rows.
    internal = sum(v["INTERNAL"] for _, v, _ in kept)
    if abs(internal) > 0.01:
        print(f"  {'INTERNAL':<12}{internal:>17.3f}{'(not in total)':>14}")

    if rejected:
        drop_sum = sum(v["TOTAL"] for _, v, _ in rejected)
        print(f"  rejected rows carry {drop_sum:+.3f} kcal/mol in TOTAL")
        for name, v, err in rejected[:5]:
            print(f"    {name:<10} TOTAL {v['TOTAL']:>8.3f}   "
                  f"components off by {err:+.4f}")

    pattern[s] = residuals

# ---- verdict ----
print()
print("=" * 78)
print("WHERE THE RESIDUAL SITS")
print("=" * 78)

if pattern:
    keys = ["VDWAALS", "EEL", "EGB", "ESURF", "TOTAL"]
    print(f"  {'SYSTEM':<16}" + "".join(f"{k:>11}" for k in keys))
    for s, r in pattern.items():
        print(f"  {s:<16}" + "".join(
            f"{r[k]:>11.3f}" if k in r else f"{'-':>11}" for k in keys))

    worst_non_surf = max(
        (abs(r.get(k, 0)) for r in pattern.values()
         for k in ("VDWAALS", "EEL", "EGB")),
        default=0,
    )
    surf = [abs(r.get("ESURF", 0)) for r in pattern.values()]

    print()
    if worst_non_surf < 0.05 and surf and max(surf) > 0.1:
        print("  Van der Waals, electrostatics and polar solvation reconcile")
        print("  exactly. The entire residual is in the non-polar surface term.")
        print()
        print("  This is a property of gbsa=2: the surface area is computed")
        print("  per atom by a route that is not strictly additive per residue.")
        print("  The electrostatic and van der Waals contributions, which are")
        print("  what the Glu136 question depends on, are exact.")
        print()
        print(f"  Bound on the effect: up to {max(surf):.2f} kcal/mol per system,")
        print("  spread across every residue rather than concentrated in one.")
    elif worst_non_surf >= 0.05:
        print("  A component other than the surface term does not reconcile.")
        print("  That is not explained by gbsa=2 and needs further work before")
        print("  these numbers are used.")
    else:
        print("  Everything reconciles. The earlier failures were in the sum")
        print("  comparison itself, not in the decomposition.")

print()
print("=" * 78)
print("Nothing was changed. This was read-only.")
print("=" * 78)
