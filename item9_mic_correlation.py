#!/usr/bin/env python3
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# ITEM 9: EXPERIMENTAL MIC vs COMPUTATIONAL PREDICTIONS
#
# Tests every prediction this study produced against the
# measured antibacterial potency:
#
#   Vina score            docking, no electrostatic term
#   MM/GBSA total         implicit solvent, charges included
#   MM/GBSA electrostatic the term that separated the compounds
#   Glu136 / FMN / His143 per-residue contributions
#   contact counts        the distance-only analysis
#
# Correlation is computed against log10(MIC), because MIC is a
# concentration and spans orders of magnitude. Lower MIC means
# more potent, and a more negative energy means tighter
# predicted binding, so a WORKING predictor gives a POSITIVE r.
#
# With three compounds no correlation is statistically
# meaningful. What can be judged is whether the rank order
# matches, so both are reported and the limitation is stated.
#
# Input: mic_values.txt in the working folder
#
#   # compound   MIC_ug_per_mL   [organism / strain]
#   5sa23        2.0             S. aureus ATCC 29213
#   17sa23       16.0            S. aureus ATCC 29213
#   20sa23       0.5             S. aureus ATCC 29213
#
#   Use >64 for "no activity at the highest concentration
#   tested"; it is treated as that value and flagged.
#
# Usage:  python3 item9_mic_correlation.py [--run focused]
# ============================================================

from __future__ import annotations

import argparse
import math
import re
from datetime import datetime
from pathlib import Path

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]
FOOTER = ("siam chowdhury | Computational and Medicinal Chemistry | "
          "Dr. Alam's Research Team (www.alamresearch.org) | Arkansas State University")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run", default="focused")
    p.add_argument("--mic", default="mic_values.txt")
    p.add_argument("--state", default="HIE",
                   help="His143 protonation state to report (HIE or HIP)")
    return p.parse_args()


args = parse_args()
run = args.run
prefix = "docking_" if run == "global" else f"docking_{run}_"
txt_dir = Path("txt_outputs") if run == "global" else Path("txt_outputs") / run
checks = Path(f"mmgbsa_checks_{run}")
step4 = Path(f"mmgbsa_{run}")

# ------------------------------------------------------------
# MIC values
# ------------------------------------------------------------
mic_path = Path(args.mic)
if not mic_path.is_file():
    raise SystemExit(
        f"ERROR: {mic_path} not found.\n\n"
        "Create it with one line per compound:\n"
        "  # compound   MIC_ug_per_mL   organism\n"
        "  5sa23        2.0             S. aureus ATCC 29213\n"
        "  17sa23       16.0            S. aureus ATCC 29213\n"
        "  20sa23       0.5             S. aureus ATCC 29213\n\n"
        "Use >64 for no activity at the highest concentration tested."
    )

mic, censored, organism = {}, set(), {}
for line in mic_path.read_text(errors="replace").splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    f = line.split()
    if len(f) < 2:
        continue
    name, value = f[0], f[1]
    if value.startswith(">"):
        censored.add(name)
        value = value[1:]
    try:
        mic[name] = float(value)
    except ValueError:
        continue
    if len(f) > 2:
        organism[name] = " ".join(f[2:])

if not mic:
    raise SystemExit(f"ERROR: no usable MIC values read from {mic_path}")

# ------------------------------------------------------------
# Computational results
# ------------------------------------------------------------
KEYS = ("VDWAALS", "EEL", "EGB", "ESURF", "TOTAL")


def read_mmgbsa(path: Path):
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


HIS_FORMS = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS"}
renum = {}
rf = step4 / "protein_amber_renum.txt"
if rf.is_file():
    for line in rf.read_text(errors="replace").splitlines():
        f = line.split()
        if len(f) >= 5 and f[2].lstrip("-").isdigit() and f[4].lstrip("-").isdigit():
            renum[int(f[4])] = f[2]


def read_residues(path: Path, wanted=("GLU136", "HIS143", "FMN")):
    """Per-residue EEL for the named residues."""
    if not path.is_file():
        return {}
    COLS = {"VDWAALS": 5, "EEL": 8, "TOTAL": 17}
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        f = line.split(",")
        if len(f) <= COLS["TOTAL"]:
            continue
        m = re.match(r"^\s*([A-Z][A-Z0-9]{1,3})\s+(\d+)\s*$", f[0])
        if not m:
            continue
        name, idx = m.group(1), int(m.group(2))
        base = HIS_FORMS.get(name, name)
        tag = "FMN" if name == "FMN" else f"{base}{renum.get(idx, idx)}"
        if tag not in wanted:
            continue
        try:
            out[tag] = {k: float(f[i]) for k, i in COLS.items()}
        except ValueError:
            pass
    return out


data = {}
for name in mic:
    row = {"mic": mic[name], "censored": name in censored}

    tsv = Path(f"{prefix}{name}/{name}_summary.tsv")
    if tsv.is_file():
        v = tsv.read_text().split("\t")
        row["vina"] = float(v[1])
        row["heavy"] = int(v[2])
        row["contacts"] = int(v[6])
        row["shared"] = int(v[7])
        row["fmn_contact"] = v[9]
        row["polar"] = int(v[10].strip())

    g = read_mmgbsa(checks / f"{name}_{args.state}" / "decomp_totals.dat") \
        or read_mmgbsa(checks / f"{name}_{args.state}" / "mmgbsa.dat") \
        or read_mmgbsa(step4 / name / "mmgbsa.dat")
    if g:
        row["gbsa"] = g

    res = read_residues(checks / f"{name}_{args.state}" / "decomp.dat")
    if res:
        row["res"] = res

    data[name] = row

usable = [n for n, r in data.items() if "vina" in r or "gbsa" in r]
if not usable:
    raise SystemExit(
        "ERROR: no computational results found for the compounds in "
        f"{mic_path}. Check that the names match the docking_* folders."
    )

# ------------------------------------------------------------
# Correlations
# ------------------------------------------------------------


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(ranks(xs), ranks(ys))


def predictor(name, getter):
    pairs = [(math.log10(data[n]["mic"]), getter(data[n]))
             for n in usable if getter(data[n]) is not None]
    if len(pairs) < 2:
        return None
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    return {
        "name": name, "n": len(pairs),
        "pearson": pearson(xs, ys), "spearman": spearman(xs, ys),
    }


def g(row, key):
    return row.get("gbsa", {}).get(key)


def r(row, tag, term="EEL"):
    return row.get("res", {}).get(tag, {}).get(term)


predictors = [
    predictor("Vina score", lambda d: d.get("vina")),
    predictor("MM/GBSA total", lambda d: g(d, "TOTAL")),
    predictor("MM/GBSA electrostatic", lambda d: g(d, "EEL")),
    predictor("MM/GBSA van der Waals", lambda d: g(d, "VDWAALS")),
    predictor("Glu136 electrostatic", lambda d: r(d, "GLU136")),
    predictor("FMN electrostatic", lambda d: r(d, "FMN")),
    predictor("His143 electrostatic", lambda d: r(d, "HIS143")),
    predictor("polar contact count", lambda d: -d["polar"] if "polar" in d else None),
]
predictors = [p for p in predictors if p]

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
NOW = f"{datetime.now():%Y-%m-%d %H:%M:%S}"
L = ["=" * 78, *IDENTITY, f"Date: {NOW}", "=" * 78, "",
     "42 EXPERIMENTAL MIC vs COMPUTATIONAL PREDICTIONS",
     "-" * 78,
     f"MIC source: {mic_path}",
     f"Docking run: {run}    His143 state used: {args.state}",
     ""]

orgs = {organism[n] for n in organism}
if orgs:
    L.append("Organism: " + "; ".join(sorted(orgs)))
    L.append("")

L += ["MEASURED AND PREDICTED",
      "-" * 78,
      f"{'COMPOUND':<11}{'MIC':>9}{'VINA':>9}{'MMGBSA':>9}{'EEL':>9}"
      f"{'GLU136':>9}{'FMN':>9}{'POLAR':>7}"]

for n in sorted(usable, key=lambda x: data[x]["mic"]):
    d = data[n]
    mv = f"{'>' if d['censored'] else ''}{d['mic']:g}"
    L.append(
        f"{n:<11}{mv:>9}"
        + f"{d.get('vina', float('nan')):>9.2f}"
        + (f"{g(d, 'TOTAL'):>9.2f}" if g(d, "TOTAL") is not None else f"{'-':>9}")
        + (f"{g(d, 'EEL'):>9.2f}" if g(d, "EEL") is not None else f"{'-':>9}")
        + (f"{r(d, 'GLU136'):>9.2f}" if r(d, "GLU136") is not None else f"{'-':>9}")
        + (f"{r(d, 'FMN'):>9.2f}" if r(d, "FMN") is not None else f"{'-':>9}")
        + (f"{d['polar']:>7}" if "polar" in d else f"{'-':>7}")
    )

L += ["",
      "MIC in ug/mL, lower is more potent. Energies in kcal/mol, more",
      "negative is tighter predicted binding. GLU136 and FMN are per-residue",
      "electrostatic contributions, where positive means repulsion.",
      ""]

# --- rank order agreement ---
by_mic = sorted(usable, key=lambda n: data[n]["mic"])
L += ["RANK ORDER", "-" * 78,
      "  experiment (most potent first):   " + " > ".join(by_mic)]

for p in predictors:
    key = {
        "Vina score": lambda d: d.get("vina"),
        "MM/GBSA total": lambda d: g(d, "TOTAL"),
        "MM/GBSA electrostatic": lambda d: g(d, "EEL"),
        "MM/GBSA van der Waals": lambda d: g(d, "VDWAALS"),
        "Glu136 electrostatic": lambda d: r(d, "GLU136"),
        "FMN electrostatic": lambda d: r(d, "FMN"),
        "His143 electrostatic": lambda d: r(d, "HIS143"),
        "polar contact count": lambda d: -d.get("polar", 0),
    }[p["name"]]
    ranked = sorted([n for n in usable if key(data[n]) is not None],
                    key=lambda n: key(data[n]))
    match = "MATCHES" if ranked == by_mic else \
            ("reversed" if ranked == by_mic[::-1] else "differs")
    L.append(f"  {p['name']:<24}  " + " > ".join(ranked) + f"   [{match}]")

L += ["", "CORRELATION WITH log10(MIC)", "-" * 78,
      "A working predictor gives a POSITIVE r: lower MIC (more potent)",
      "together with a more negative energy.",
      "",
      f"{'PREDICTOR':<26}{'n':>4}{'PEARSON':>10}{'SPEARMAN':>11}   DIRECTION"]

for p in predictors:
    pr, sp = p["pearson"], p["spearman"]
    if math.isnan(pr):
        direction = "undefined"
    elif pr > 0.3:
        direction = "as expected"
    elif pr < -0.3:
        direction = "OPPOSITE to expected"
    else:
        direction = "no relationship"
    L.append(f"{p['name']:<26}{p['n']:>4}{pr:>10.3f}{sp:>11.3f}   {direction}")

n_c = len(usable)
L += ["",
      "STATISTICAL WEIGHT",
      "-" * 78,
      f"n = {n_c}. "
      + ("With three compounds a correlation coefficient carries almost no"
         if n_c <= 3 else "At this sample size a correlation carries limited"),
      "statistical weight: r can reach 1.0 by chance, and no p-value would be",
      "meaningful. Read the RANK ORDER section as the real test, and treat the",
      "correlations as description rather than evidence. A conclusion about",
      "which predictor works needs the full compound series from the paper.",
      ""]

if censored:
    L += ["Censored values (>) were treated as equal to the stated limit,",
          f"which understates the potency gap for: {', '.join(sorted(censored))}",
          ""]

# --- interpretation ---
L += ["=" * 78, "WHAT THIS MEANS", "=" * 78]

best = max((p for p in predictors if not math.isnan(p["pearson"])),
           key=lambda p: p["pearson"], default=None)
worst = min((p for p in predictors if not math.isnan(p["pearson"])),
            key=lambda p: p["pearson"], default=None)

if best and best["pearson"] > 0.3:
    L += [f"  Best agreement: {best['name']} (r = {best['pearson']:.2f}).",
          "  This is consistent with the compounds acting on FabK in the way",
          "  modelled here, though it does not establish it.",
          ""]
if worst and worst["pearson"] < -0.3:
    L += [f"  Opposite to expectation: {worst['name']} "
          f"(r = {worst['pearson']:.2f}).",
          "  A predictor pointing the wrong way is more informative than a",
          "  weak one: it suggests the quantity it measures is not what",
          "  drives potency in this series.",
          ""]
if best and best["pearson"] <= 0.3 and (not worst or worst["pearson"] >= -0.3):
    L += ["  No predictor tracks MIC in either direction.",
          "",
          "  The most likely explanation is that these compounds do not",
          "  inhibit FabK by binding this pocket. The MM/GBSA showed the site",
          "  is strongly anionic (Glu136 and the FMN phosphate) while these",
          "  compounds are benzoate anions, so all three are electrostatically",
          "  repelled, and the crystal inhibitor XCJ scores far better than any",
          "  of them. The paper assigned the fatty-acid pathway by CRISPRi,",
          "  which identifies the pathway rather than the enzyme.",
          ""]

L += ["  Regardless of the outcome above, the docking protocol itself was",
      "  validated: XCJ redocked to 0.82 A RMSD and scored favourably on",
      "  every measure. A protocol can be correct and still show that a",
      "  compound series does not fit the target it was tested against.",
      ""]

L += ["=" * 78, "LIMITATIONS", "=" * 78,
      "  MIC measures whole-cell antibacterial activity. It depends on",
      "  membrane permeability, efflux, metabolic stability and off-target",
      "  effects as well as target binding, so even a correct binding model",
      "  need not track MIC.",
      "  Single-structure MM/GBSA, no entropy, one docked pose per compound.",
      "  Vina and MM/GBSA values are computational scores, not experimentally",
      "  measured binding free energies.",
      ""]

txt_dir.mkdir(parents=True, exist_ok=True)
out = txt_dir / "42_MIC_vs_predictions.txt"
out.write_text("\n".join(L) + "\n")
print("\n".join(L[5:]))
print(f"\nWritten: {out}")

# ------------------------------------------------------------
# Graphs
# ------------------------------------------------------------
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


panels = [
    ("Vina score (kcal/mol)", lambda d: d.get("vina")),
    ("MM/GBSA total (kcal/mol)", lambda d: g(d, "TOTAL")),
    ("MM/GBSA electrostatic (kcal/mol)", lambda d: g(d, "EEL")),
    ("Glu136 electrostatic (kcal/mol)", lambda d: r(d, "GLU136")),
]
panels = [(t, f) for t, f in panels
          if sum(1 for n in usable if f(data[n]) is not None) >= 2]

if panels:
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5))
    if len(panels) == 1:
        axes = [axes]
    for ax, (title, f) in zip(axes, panels):
        for n in usable:
            v = f(data[n])
            if v is None:
                continue
            ax.scatter(data[n]["mic"], v, s=80,
                       marker="v" if data[n]["censored"] else "o")
            ax.annotate(n, (data[n]["mic"], v), xytext=(6, 4),
                        textcoords="offset points", fontsize=9)
        ax.set_xscale("log")
        ax.set_xlabel("MIC (ug/mL, log scale)")
        ax.set_ylabel(title)
        ax.grid(True, linestyle="--", alpha=0.4)
    fig.suptitle(f"Experimental potency vs predictions: 7L00 FabK [{run}]")
    finish(fig, graph_dir / "MIC_vs_all_predictions.png")
    print(f"Graph: {graph_dir}/MIC_vs_all_predictions.png")

if len(predictors) >= 2:
    ps = [p for p in predictors if not math.isnan(p["pearson"])]
    ps.sort(key=lambda p: p["pearson"])
    fig = plt.figure(figsize=(9, max(4, 0.55 * len(ps) + 2)))
    colors = ["#C44E52" if p["pearson"] < 0 else "#55A868" for p in ps]
    plt.barh([p["name"] for p in ps], [p["pearson"] for p in ps], color=colors)
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Pearson r against log10(MIC)   (positive = as expected)")
    plt.title(f"Which prediction tracks potency?  (n = {n_c})")
    plt.xlim(-1.05, 1.05)
    plt.grid(True, axis="x", linestyle="--", alpha=0.4)
    finish(fig, graph_dir / "MIC_predictor_comparison.png")
    print(f"Graph: {graph_dir}/MIC_predictor_comparison.png")
