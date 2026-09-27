#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# EXPORT WORKFLOW
#
# Writes the entire project into one plain-text file: the
# commands in the order they were run, the full source of every
# script, the inputs, and every report produced. The result is
# a complete record of the study that can be re-read, reused as
# a template for the next target, or handed to someone else as
# a worked example.
#
# Nothing is moved, renamed or deleted.
#
#   ./export_workflow.sh                  full export
#   ./export_workflow.sh --code-only      scripts and commands, no output
#   ./export_workflow.sh -o myfile.txt    choose the file name
# ============================================================

set -Eeuo pipefail

OUT="7L00_FabK_complete_workflow.txt"
CODE_ONLY=0

while (( $# )); do
    case "$1" in
        --code-only) CODE_ONLY=1; shift ;;
        -o) OUT="$2"; shift 2 ;;
        *) echo "unknown option: $1"; exit 1 ;;
    esac
done

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found."; exit 1; }

python3 - "$OUT" "$CODE_ONLY" <<'PY'
import sys
from datetime import datetime
from pathlib import Path

out_path = Path(sys.argv[1])
code_only = sys.argv[2] == "1"

IDENTITY = [
    "Siam Chowdhury (www.siamchowdhury.com)",
    "Computational and Medicinal Chemistry",
    "Mentor: Mohammad Alam PhD., Professor of Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University, AR, USA",
]

W = 78
BAR = "=" * W
THIN = "-" * W


def head(title, level=1):
    if level == 1:
        return ["", BAR, BAR, f"  {title}", BAR, BAR, ""]
    if level == 2:
        return ["", BAR, f"  {title}", BAR, ""]
    return ["", THIN, f"  {title}", THIN, ""]


def block(path: Path, label=""):
    """A file's full contents, fenced so the boundaries are obvious."""
    text = path.read_text(errors="replace").rstrip()
    n = len(text.splitlines())
    return [
        f"<<<<<< BEGIN {label or path.name}   ({n} lines, {path.stat().st_size:,} bytes)",
        "",
        text,
        "",
        f">>>>>> END {label or path.name}",
        "",
    ]


# ============================================================
# The workflow, stage by stage
# ============================================================
STAGES = [
    {
        "title": "STAGE 1 - STRUCTURE PREPARATION",
        "what": [
            "Download PDB 7L00 (C. difficile FabK with inhibitor XCJ and",
            "cofactor FMN) and split it into protein, FMN, XCJ, water and",
            "other components. No graphical program is used at any point:",
            "everything runs from bash and Python so the work is reproducible.",
        ],
        "scripts": ["prepare_7L00_structure.sh"],
        "commands": [
            "cd ~/chem/s_series",
            "conda activate docking",
            "mkdir -p txt_outputs",
            "chmod +x prepare_7L00_structure.sh",
            "./prepare_7L00_structure.sh 2>&1 | tee txt_outputs/00_step1_terminal.txt",
        ],
        "reports": ["01_", "02_", "03_", "04_", "05_", "06_"],
        "notes": [
            "7L00 has four chains, A to D, each with one XCJ and one FMN.",
            "Report 06 pairs each XCJ with its nearest FMN, which is how the",
            "pocket for docking was chosen.",
        ],
    },
    {
        "title": "STAGE 2 - RECEPTOR, GRID AND REDOCKING VALIDATION",
        "what": [
            "Build the receptor from chains A and B plus both FMN molecules,",
            "characterise the binding site, set the grid on the crystal XCJ,",
            "then redock XCJ and measure the RMSD against its crystal pose.",
            "This is the control that decides whether the protocol can be",
            "trusted on unknown compounds.",
        ],
        "scripts": ["step2_7L00_receptor_redock.sh"],
        "commands": [
            "chmod +x step2_7L00_receptor_redock.sh",
            "./step2_7L00_receptor_redock.sh 2>&1 | tee txt_outputs/00_step2_terminal.txt",
        ],
        "reports": ["07_", "08_", "09_", "10_", "11_", "12_", "13_", "14_",
                    "16_", "17_", "17b_", "18_"],
        "notes": [
            "RESULT: RMSD 0.82 A, well inside the 2 A threshold. XCJ scored",
            "-10.11 kcal/mol, the reference for every later comparison.",
            "",
            "Two details that mattered:",
            "  - Ligand bond orders come from the RCSB chemical dictionary,",
            "    not guessed from a PDB file, which has no bond orders.",
            "  - Residues with incomplete side chains are trimmed to ALA or",
            "    GLY rather than deleted, because Meeko rejects partial",
            "    residues and would otherwise drop them entirely.",
        ],
    },
    {
        "title": "STAGE 3 - DOCKING THE COMPOUNDS",
        "what": [
            "Prepare each compound at pH 7.4, dock it with the validated",
            "protocol, export poses and complexes, analyse contacts, rank and",
            "graph. Run once with a wide box and once with a narrow box.",
        ],
        "scripts": ["step3_7L00_dock_compounds.sh", "step3b_pose_scan.py"],
        "commands": [
            "chmod +x step3_7L00_dock_compounds.sh",
            "",
            "# wide box (22 A), all *.mol files in the folder",
            "./step3_7L00_dock_compounds.sh 2>&1 | tee txt_outputs/00_step3_terminal.txt",
            "",
            "# narrow box (16 A), forced into the FMN pocket",
            "RUN=focused BOX=16 ./step3_7L00_dock_compounds.sh \\",
            "    2>&1 | tee txt_outputs/00_step3_focused_terminal.txt",
            "",
            "# redraw graphs only, no redocking",
            "GRAPHS_ONLY=1 ./step3_7L00_dock_compounds.sh",
            "",
            "# anchor-residue contacts across all 20 poses, not just the best",
            "python3 step3b_pose_scan.py --run focused",
            "python3 step3b_pose_scan.py",
        ],
        "reports": ["19_", "20_", "21_", "23_", "24_", "25_", "26_", "27_", "34_"],
        "notes": [
            "In the wide box all three compounds settled ~9.7 A from the XCJ",
            "position and never touched FMN. The narrow box forced them into",
            "the real pocket at a cost of 1.1 to 2.1 kcal/mol. The focused",
            "poses were used for everything afterwards.",
            "",
            "The three docking scores fell within 0.8 kcal/mol of each other,",
            "which is inside Vina's own uncertainty, so docking alone could",
            "not rank them. That is why the study continued to MM/GBSA.",
        ],
    },
    {
        "title": "STAGE 4 - MM/GBSA RESCORING",
        "what": [
            "Assign GAFF2 parameters and AM1-BCC charges to each ligand and to",
            "FMN, build complex, receptor and ligand topologies with ff14SB,",
            "minimise with the backbone restrained, and run MM/GBSA.",
        ],
        "scripts": ["step4_mmgbsa.sh"],
        "commands": [
            "# install AmberTools once, if not already present",
            "# conda install -c conda-forge ambertools -y",
            "",
            "chmod +x step4_mmgbsa.sh",
            "MAXCYC=500 NCYC=200 ./step4_mmgbsa.sh \\",
            "    2>&1 | tee txt_outputs/focused/00_step4_terminal.txt",
        ],
        "reports": ["35_", "36_", "37_", "38_"],
        "notes": [
            "Three failures had to be fixed here, each worth remembering:",
            "",
            "  1. antechamber's bondtype crashed on the FMN phosphate. Cause:",
            "     it was given a PDB file, which carries no bond orders, so it",
            "     had to guess them. Fix: feed it the RCSB SDF instead.",
            "",
            "  2. sqm refused to run with an odd electron count. Cause: the",
            "     ligand charge was read as 0 when these carboxylates are -1.",
            "     Fix: read formal charges with RDKit and check that the",
            "     electron count is even BEFORE calling antechamber.",
            "",
            "  3. A MOL2 conversion silently dropped the formal charge while",
            "     the SDF kept it. Fix: trust the SDF, and let electron parity",
            "     settle any disagreement between sources.",
            "",
            "RESULT: 20sa23 -32.02, 5sa23 -17.89, 17sa23 -10.78 kcal/mol,",
            "against -43.69 for XCJ. Every EEL term for the compounds was",
            "POSITIVE, meaning electrostatics opposed binding. That unexpected",
            "result is what forced stage 5.",
        ],
    },
    {
        "title": "STAGE 5 - VALIDATING THE MM/GBSA RESULT",
        "what": [
            "The MM/GBSA ranking disagreed with the docking scores and every",
            "electrostatic term came out positive. Before reporting it, three",
            "checks tested whether the result was real or an artefact:",
            "",
            "  A  control      MM/GBSA on the redocked XCJ pose, known correct",
            "                  to 0.82 A. If a correct pose also gave positive",
            "                  electrostatics, the fault would be the setup.",
            "  B  protonation  His143 is catalytic and the preparation assumed",
            "                  the neutral form. Both states were tested.",
            "  C  decomposition  per-residue energies, with a consistency",
            "                  check on the decomposition itself.",
        ],
        "scripts": ["step6_mmgbsa_checks.sh", "step6c_decomposition.sh",
                    "diagnose_decomp.sh", "diagnose_decomp_residual.py"],
        "commands": [
            "chmod +x step6_mmgbsa_checks.sh step6c_decomposition.sh",
            "",
            "# XCJ control plus both His143 states: 8 systems",
            "./step6_mmgbsa_checks.sh 2>&1 | tee txt_outputs/focused/00_step6_terminal.txt",
            "",
            "# per-residue decomposition, verified",
            "./step6c_decomposition.sh 2>&1 | tee txt_outputs/focused/00_step6c_terminal.txt",
            "",
            "# read-only diagnostics used when decomposition first failed",
            "./diagnose_decomp.sh",
            "python3 diagnose_decomp_residual.py",
        ],
        "reports": ["39_", "40_", "41_"],
        "notes": [
            "A  XCJ gives EEL = -42.51, favourable. The positive values for",
            "   the compounds are therefore real properties of the compounds.",
            "",
            "B  Energies shift 2 to 6 kcal/mol between HIE and HIP, and the",
            "   compound order does not change. The ranking does not depend on",
            "   that assumption.",
            "",
            "C  Two failures here, both instructive:",
            "     - 'Mismatch in number of decomp terms': print_res listed the",
            "       protein residues but not the ligand, so the complex and",
            "       ligand decompositions could never match. Fix: request",
            "       every residue.",
            "     - Sums disagreed with totals by 20 to 37 kcal/mol. Cause was",
            "       the PARSER, not AmberTools: a decimal regex matched",
            "       '8.88e-16' as 8.88, discarding the exponent. Fix: parse the",
            "       CSV by column position.",
            "",
            "   After the fix, VDW, EEL and EGB reconcile to within 0.014",
            "   kcal/mol. ESURF does not decompose exactly, because gbsa=2",
            "   computes surface area per atom by a route that is not additive",
            "   per residue; that residual is bounded at 1.10 kcal/mol and",
            "   reported rather than hidden.",
            "",
            "RESULT: the pocket is strongly anionic. Glu136 (-1) and the FMN",
            "phosphate (-2) repel the compounds' carboxylate: Glu136 gives",
            "+12.7 to +34.1 kcal/mol of repulsion and FMN +28.3 to +49.8,",
            "while for XCJ Glu136 is attractive at -11.0. The one favourable",
            "anchor is His143, -19.6 to -26.3 when protonated, not enough to",
            "overcome the rest.",
        ],
    },
    {
        "title": "STAGE 6 - COMPARISON WITH MEASURED MIC",
        "what": [
            "Test every prediction the study produced against the measured",
            "antibacterial potency reported in the source paper.",
        ],
        "scripts": ["item9_mic_correlation.py"],
        "commands": [
            "# create mic_values.txt first:",
            "#   # compound   MIC_ug_per_mL   organism",
            "#   5sa23        16.0            S. aureus ATCC 29213",
            "#   17sa23       32.0            S. aureus ATCC 29213",
            "#   20sa23       1.0             S. aureus ATCC 29213",
            "",
            "python3 item9_mic_correlation.py \\",
            "    2>&1 | tee txt_outputs/focused/00_item9_terminal.txt",
        ],
        "reports": ["42_"],
        "notes": [
            "RESULT: the MM/GBSA order matches the measured order exactly:",
            "20sa23 (MIC 1) > 5sa23 (16) > 17sa23 (32).",
            "",
            "The docking score did not match. The count of close polar",
            "contacts gave the order BACKWARDS, because counting contacts by",
            "distance alone treats a carboxylate pressed against Glu136 as a",
            "hydrogen bond when both carry a negative charge and repel. The",
            "Vina scoring function has no electrostatic term either, which is",
            "why it gave no useful signal.",
            "",
            "CAUTION: with three compounds a correct rank order has a",
            "one-in-six chance of arising by luck. What supports it is the",
            "mechanism, which is testable by docking more compounds.",
        ],
    },
    {
        "title": "STAGE 7 - REPORTING AND PUBLISHING",
        "what": [
            "Assemble everything into one self-contained HTML page and push",
            "the project to GitHub with a live page.",
        ],
        "scripts": ["build_report.sh", "show_all_results.sh", "export_workflow.sh"],
        "commands": [
            "chmod +x build_report.sh",
            "./build_report.sh",
            "",
            "# publish",
            "cp 7L00_FabK_docking_study.html index.html",
            "git add -A",
            "git commit -m 'describe what changed'",
            "git push",
            "",
            "# this file",
            "./export_workflow.sh",
        ],
        "reports": [],
        "notes": [
            "Repository: https://github.com/siamchowdhury-1/fabk-docking-7L00",
            "Live report: https://siamchowdhury-1.github.io/fabk-docking-7L00/",
        ],
    },
]

# ============================================================
# Assemble
# ============================================================
root = Path(".")
L = []

L += [
    BAR, BAR,
    "",
    "  COMPLETE WORKFLOW RECORD",
    "  Molecular docking and MM/GBSA study",
    "  Pyrazole benzoic acid derivatives vs C. difficile FabK (PDB 7L00)",
    "",
    BAR,
    "",
]
L += [f"  {line}" for line in IDENTITY]
L += [
    "",
    f"  Exported: {datetime.now():%d %B %Y, %H:%M}",
    f"  Working directory: {root.resolve()}",
    "",
    BAR, BAR,
]

L += head("WHAT THIS FILE IS", 2)
L += [
    "A complete record of one computational study: the commands in the order",
    "they were run, the full source of every script, and every report the",
    "scripts produced. It can be read as a worked example, reused as a",
    "template for a different target, or handed to someone else so they can",
    "repeat the work.",
    "",
    "To adapt it to a new target, the things that change are:",
    "  - the PDB ID and the name of the co-crystallised ligand",
    "  - the cofactor, if there is one, and its formal charge",
    "  - the chain and residue number of the reference ligand",
    "  - the ligand .mol files",
    "  - the catalytic residue tested for protonation",
    "Everything else is general.",
    "",
]

L += head("HEADLINE RESULT", 2)
L += [
    "Computational study of fatty acid biosynthesis inhibition at FabK",
    "(PDB 7L00). MM/GBSA reproduced the measured antibacterial potency order",
    "exactly:",
    "",
    "    20sa23 (MIC 1 ug/mL)  >  5sa23 (MIC 16 ug/mL)  >  17sa23 (MIC 32 ug/mL)",
    "",
    "matching the MIC values reported in the source paper. The docking score",
    "alone did not, and the count of close polar contacts gave the order",
    "backwards.",
    "",
    "Source paper: Development and Antibacterial Properties of",
    "4-[4-(Anilinomethyl)-3-phenylpyrazol-1-yl]benzoic Acid Derivatives as",
    "Fatty Acid Biosynthesis Inhibitors",
    "https://doi.org/10.1021/acs.jmedchem.3c00969",
    "",
]

L += head("SOFTWARE USED", 2)
L += [
    "  conda environment : docking",
    "  AutoDock Vina     : docking",
    "  Meeko             : receptor and ligand PDBQT preparation",
    "  Open Babel        : 3D generation, protonation, format conversion",
    "  RDKit             : formal charges, symmetry-aware RMSD",
    "  AmberTools        : antechamber, parmchk2, tleap, sander, MMPBSA.py,",
    "                      cpptraj, pdb4amber",
    "  matplotlib        : graphs",
    "",
    "  Everything runs from the command line. No graphical program is used",
    "  at any stage, so every step is scripted and reproducible.",
    "",
]

L += head("CONTENTS", 2)
for i, st in enumerate(STAGES, start=1):
    L.append(f"  {i}. {st['title']}")
L += [
    f"  {len(STAGES) + 1}. INPUT FILES",
    f"  {len(STAGES) + 2}. FULL SCRIPT SOURCE",
]
if not code_only:
    L.append(f"  {len(STAGES) + 3}. ALL OUTPUT REPORTS")
L.append("")

# ---- stages ----
txt_root = Path("txt_outputs")
used_reports = set()

for i, st in enumerate(STAGES, start=1):
    L += head(f"{i}. {st['title']}", 1)
    L += ["PURPOSE", THIN] + st["what"] + [""]

    if st["scripts"]:
        L += ["SCRIPTS USED", THIN]
        for s in st["scripts"]:
            mark = "" if Path(s).is_file() else "   (not present in this folder)"
            L.append(f"  {s}{mark}")
        L.append("")

    L += ["COMMANDS", THIN] + [f"  {c}" if c else "" for c in st["commands"]] + [""]

    if st["notes"]:
        L += ["NOTES AND RESULTS", THIN] + st["notes"] + [""]

    if st["reports"] and not code_only:
        files = sorted(
            p for p in txt_root.rglob("*.txt")
            if p.is_file() and p.stat().st_size > 0
            and any(p.name.startswith(pre) for pre in st["reports"])
        )
        if files:
            L += ["REPORTS PRODUCED", THIN]
            for p in files:
                L.append(f"  {p}")
                used_reports.add(p)
            L.append("")

# ---- inputs ----
n = len(STAGES) + 1
L += head(f"{n}. INPUT FILES", 1)

mols = sorted(root.glob("*.mol"))
if mols:
    L += ["Compound structures as drawn:", ""]
    for p in mols:
        L += block(p)

for extra in ("mic_values.txt", "7L00_AB_grid_box.txt"):
    p = Path(extra)
    if p.is_file():
        L += [f"{extra}:", ""] + block(p)

L += [
    "The receptor structure 7L00.pdb is downloaded by stage 1 and is not",
    "reproduced here because of its size. It is available from the Protein",
    "Data Bank at https://files.rcsb.org/download/7L00.pdb",
    "",
]

# ---- scripts ----
n += 1
L += head(f"{n}. FULL SCRIPT SOURCE", 1)
L += [
    "Every script in full, in the order it is used. Each is self-contained",
    "and documents its own options at the top.",
    "",
]

ordered, seen = [], set()
for st in STAGES:
    for s in st["scripts"]:
        if s not in seen and Path(s).is_file():
            ordered.append(s)
            seen.add(s)
for p in sorted(root.glob("*.sh")) + sorted(root.glob("*.py")):
    if p.name not in seen:
        ordered.append(p.name)
        seen.add(p.name)

for s in ordered:
    L += head(s, 3)
    L += block(Path(s))

# ---- output ----
if not code_only:
    n += 1
    L += head(f"{n}. ALL OUTPUT REPORTS", 1)
    L += [
        "Every report the scripts produced, in full. Each records what the",
        "stage did, what it produced and any warnings, so that every number",
        "in the study can be traced back to the run that produced it.",
        "",
    ]

    all_reports = sorted(
        p for p in txt_root.rglob("*.txt")
        if p.is_file() and p.stat().st_size > 0
        and not p.name.startswith("00_")
    )

    group = {}
    for p in all_reports:
        key = "16 A box (focused run)" if "focused" in p.parts else "22 A box (global run)"
        group.setdefault(key, []).append(p)

    for key in sorted(group):
        L += head(key, 2)
        for p in group[key]:
            L += head(str(p), 3)
            L += block(p, str(p))

L += [
    "", BAR, BAR, "",
    "  END OF WORKFLOW RECORD",
    "",
]
L += [f"  {line}" for line in IDENTITY]
L += [
    "",
    "  Docking and MM/GBSA values are computational scores, not",
    "  experimentally measured binding free energies.",
    "",
    BAR, BAR,
]

out_path.write_text("\n".join(L) + "\n", encoding="utf-8")
size = out_path.stat().st_size
lines = len(L)
print(f"{lines:,} lines, {size / 1024:.0f} KB -> {out_path}")
print(f"  scripts included : {len(ordered)}")
if not code_only:
    print(f"  reports included : {len(all_reports)}")
PY

echo
echo "Written: $(pwd)/$OUT"
echo "Read with:  less $OUT"
