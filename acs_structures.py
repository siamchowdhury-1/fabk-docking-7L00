#!/usr/bin/env python3
# Siam Chowdhury | [Dr. Alam's Research Team] | Arkansas State University
# ACS 1996 style 2D structures for the five-compound series.
# Reads only the .mol files. No score, report or structure is touched.
import sys
from pathlib import Path
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdFMCS, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:
    sys.exit("ERROR: RDKit not found. Run: conda activate docking")
RDLogger.DisableLog("rdApp.*")

SERIES = [
    ("20sa23", "3,5-dichloro",  1, -32.02),
    ("25sa23", "3-CF3, 5-F",    2, -28.53),
    ("22sa23", "3-CF3, 4-F",    4, -23.35),
    ("5sa23",  "3-fluoro",     16, -17.89),
    ("17sa23", "3,4-difluoro", 32, -10.78),
]
args = sys.argv[1:]
ALIGN = "--no-align" not in args
SCALE = int(args[args.index("--scale") + 1]) if "--scale" in args else 3
OUT = Path(args[args.index("--out") + 1]).expanduser() if "--out" in args \
      else Path("acs_structures")
BOND_PT = 14.4
W, H = 340 * SCALE, 250 * SCALE

print("=" * 62)
print("siam chowdhury")
print("[Dr. Alam's Research Team] www.alamresearch.org")
print("Arkansas State University")
print("=" * 62)
print()
print("ACS-STYLE 2D STRUCTURES")
print("Alignment on common core: %s" % ("yes" if ALIGN else "no"))
print("Canvas: %d x %d px   bond length: %.1f pt" % (W, H, BOND_PT * SCALE))
print()

mols, missing = [], []
for name, aniline, mic, gbsa in SERIES:
    p = Path(name + ".mol")
    if not p.is_file():
        missing.append(name); continue
    m = Chem.MolFromMolFile(str(p))
    if m is None:
        missing.append(name + " (unreadable)"); continue
    m.SetProp("_Name", name)
    mols.append((name, aniline, mic, gbsa, m))
if missing:
    print("MISSING OR UNREADABLE")
    for n in missing: print("  " + n)
    print()
if not mols:
    sys.exit("ERROR: no structures read. Run from the study directory.")

if ALIGN and len(mols) > 1:
    raw = [m for *_, m in mols]
    mcs = rdFMCS.FindMCS(raw, timeout=20, ringMatchesRingOnly=True,
                         completeRingsOnly=True)
    core = Chem.MolFromSmarts(mcs.smartsString)
    if core is not None and core.GetNumAtoms() > 5:
        AllChem.Compute2DCoords(core)
        for *_, m in mols:
            AllChem.GenerateDepictionMatching2DStructure(m, core,
                                                         acceptFailure=True)
        print("Common core: %d atoms, all structures aligned to it"
              % core.GetNumAtoms())
    else:
        for *_, m in mols: AllChem.Compute2DCoords(m)
        print("No usable common core found; laid out independently")
else:
    for *_, m in mols: AllChem.Compute2DCoords(m)
print()

OUT.mkdir(parents=True, exist_ok=True)

def draw(mol, stem, legend):
    prepared = rdMolDraw2D.PrepareMolForDrawing(mol)
    mean_bond = rdMolDraw2D.MeanBondLength(prepared)
    for cls, suffix, binary in ((rdMolDraw2D.MolDraw2DSVG, ".svg", False),
                                (rdMolDraw2D.MolDraw2DCairo, ".png", True)):
        d = cls(W, H)
        o = d.drawOptions()
        rdMolDraw2D.SetACS1996Mode(o, mean_bond * SCALE)
        o.fixedBondLength = BOND_PT * SCALE
        if legend: o.legendFontSize = int(11 * SCALE)
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, legend=legend)
        d.FinishDrawing()
        t = d.GetDrawingText()
        Path(str(stem) + suffix).write_bytes(t if binary else t.encode())

print("WRITING")
print("-" * 62)
for name, aniline, mic, gbsa, m in mols:
    legend = "%s   %s   MIC %s ug/mL   %.2f kcal/mol" % (name, aniline, mic, gbsa)
    draw(m, OUT / name, legend)
    print("  %-10s %s.svg  +  .png   (%d heavy atoms)"
          % (name, OUT / name, m.GetNumAtoms()))

if len(mols) > 1:
    legends = ["%s  (%s)  MIC %s" % (n, a, mic) for n, a, mic, _, _ in mols]
    img = Draw.MolsToGridImage([m for *_, m in mols], molsPerRow=3,
                               subImgSize=(int(330 * SCALE / 2),
                                           int(260 * SCALE / 2)),
                               legends=legends, useSVG=True)
    svg = img.data if hasattr(img, "data") else img
    (OUT / "series_panel.svg").write_text(svg)
    print("  %-10s %s   all %d compounds, ordered by potency"
          % ("panel", OUT / "series_panel.svg", len(mols)))

print()
print("-" * 62)
print("Use the SVG for the manuscript: vector, sharp at any size.")
print("ACS single-column width is 3.33 in, double-column 7 in.")
print("N-H and O-H are labelled explicitly; C-H stays implicit.")
