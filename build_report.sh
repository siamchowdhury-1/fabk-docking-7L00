#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# GitHub: https://github.com/siamchowdhury-1
# ORCID:  https://orcid.org/0000-0001-5216-6879
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# Mail:   malam@astate.edu
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# BUILD PROJECT REPORT
#
# Assembles one self-contained HTML page holding the whole
# study: methods with definitions, every script, every text
# report, every graph, and the final comparison against the
# measured MIC values. Images and code are embedded, so the
# file stands alone and can be published as-is.
#
#   ./build_report.sh              build and open
#   ./build_report.sh              build only
# ============================================================

set -Eeuo pipefail

OUT="7L00_FabK_docking_study.html"
OPEN=1
for arg in "$@"; do
    case "$arg" in
        --no-open) OPEN=0 ;;
        *) echo "unknown option: $arg"; exit 1 ;;
    esac
done

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found."; exit 1; }

echo "Assembling report in $(pwd) ..."

python3 - "$OUT" <<'PY'
import base64
import html
import re
import sys
from datetime import datetime
from pathlib import Path

out_path = Path(sys.argv[1])

AUTHOR = "Siam Chowdhury"
AUTHOR_SITE = "www.siamchowdhury.com"
AUTHOR_GITHUB = "https://github.com/siamchowdhury-1"
AUTHOR_ORCID = "https://orcid.org/0000-0001-5216-6879"
MENTOR_MAIL = "malam@astate.edu"
AFFIL = [
    "Computational and Medicinal Chemistry",
    "Mentor: Mohammad Alam, PhD &mdash; Professor of Chemistry "
    f"(<a href=\"mailto:{MENTOR_MAIL}\">{MENTOR_MAIL}</a>)",
    "Dr. Alam's Research Team &mdash; "
    "<a href=\"https://www.alamresearch.org\">www.alamresearch.org</a>",
    "Arkansas State University, AR, USA",
]
NOW = f"{datetime.now():%d %B %Y}"

PAPER_TITLE = ("Development and Antibacterial Properties of "
               "4-[4-(Anilinomethyl)-3-phenylpyrazol-1-yl]benzoic Acid "
               "Derivatives as Fatty Acid Biosynthesis Inhibitors")
PAPER_DOI = "https://doi.org/10.1021/acs.jmedchem.3c00969"

# ============================================================
# Measured and computed results
# ============================================================
COMPOUNDS = [
    # name, aniline substitution, MIC, Vina, MM/GBSA, EEL,
    # Glu136 EEL, FMN EEL, polar contacts    (sorted by measured MIC)
    ("20sa23", "3,5-dichloro", "1",  -7.30, -32.02,  55.49, 12.73, 28.27, 2),
    ("25sa23", "3-CF3, 5-F",   "2",  -8.19, -28.53,  67.76, 14.20, 33.49, 2),
    ("22sa23", "3-CF3, 4-F",   "4",  -8.36, -23.35, 111.19, 32.46, 36.14, 6),
    ("5sa23",  "3-fluoro",     "16", -7.99, -17.89, 121.59, 34.08, 43.24, 5),
    ("17sa23", "3,4-difluoro", "32", -7.17, -10.78, 125.80, 32.70, 49.82, 8),
]
XCJ = ("XCJ (crystal inhibitor)", "co-crystallised", "-",
       -10.11, -43.69, -42.51, -11.03, 2.38, "-")

# MM/GBSA energy components, neutral His143 (kcal/mol), by measured potency.
COMPONENTS = [
    # name, VDWAALS, EEL, EGB, ESURF, TOTAL
    ("20sa23", -45.74,  55.49, -36.68, -5.09, -32.02),
    ("25sa23", -47.37,  67.76, -42.96, -5.96, -28.53),
    ("22sa23", -59.94, 111.19, -68.22, -6.38, -23.35),
    ("5sa23",  -54.46, 121.59, -79.36, -5.67, -17.89),
    ("17sa23", -51.36, 125.80, -79.66, -5.56, -10.78),
    ("XCJ",    -60.97, -42.51,  65.89, -6.10, -43.69),
]

# Per-residue electrostatic contributions, neutral His143 unless noted.
# Each entry: name, Glu136 EEL, His143 EEL (HIE), His143 EEL (HIP), FMN EEL
RESIDUES = [
    ("20sa23", 12.73,  -3.50, -26.34, 28.27),
    ("25sa23", 14.20,  -1.90, -28.58, 33.49),
    ("22sa23", 32.46,   5.04, -21.79, 36.14),
    ("5sa23",  34.08,   5.52, -24.95, 43.24),
    ("17sa23", 32.70,   7.02, -19.59, 49.82),
    ("XCJ",   -11.03,  -3.60,  -7.32,  2.38),
]

# Total binding energy in each protonation state (kcal/mol).
PROTONATION = [
    # name, HIE, HIP
    ("20sa23", -32.02, -34.13),
    ("25sa23", -28.53, -30.60),
    ("22sa23", -23.35, -29.66),
    ("5sa23",  -17.89, -21.62),
    ("17sa23", -10.78, -14.36),
    ("XCJ",    -43.69, -49.88),
]

# ============================================================
# Methods: number, title, and the blocks that make it up
# ============================================================
# Sections still to be carried out. Flagged in the page rather than hidden,
# so the report states plainly what has and has not been done.
PENDING = {7, 8}

METHODS = [
    (1, "Structure preparation", [
        ("Definitions",
         "<b>PDB</b> = Protein Data Bank. <b>PDB ID 7L00</b> = crystal structure "
         "of <i>Clostridioides difficile</i> FabK in complex with the inhibitor "
         "XCJ and the cofactor FMN, solved at 1.72 &Aring; resolution. "
         "<b>FabK</b> = enoyl-acyl carrier protein reductase II, an enzyme of "
         "bacterial fatty acid biosynthesis. <b>FMN</b> = flavin mononucleotide, "
         "the redox cofactor FabK requires. <b>XCJ</b> = the three-letter "
         "chemical component code of the co-crystallised inhibitor. "
         "<b>&Aring;</b> = &aring;ngstr&ouml;m, 10<sup>-10</sup> m. "
         "<b>HETATM</b> = the record type a PDB file uses for atoms that are "
         "not part of the standard protein or nucleic acid."),
        ("What I did",
         "I downloaded 7L00 and separated it into its components with a script "
         "rather than a graphical program, so that every step is reproducible. "
         "I kept the protein and FMN, and removed the co-crystallised inhibitor "
         "XCJ, the crystallographic waters, glycerol (GOL) and sodium ions (NA). "
         "The structure contains four chains, A to D, each with one XCJ and one "
         "FMN. I built the receptor from chains A and B, because the pocket of "
         "chain A lies close to chain B."),
        ("Why",
         "FMN is part of the active site architecture, so removing it would "
         "distort the pocket and give unrealistic poses. XCJ has to come out to "
         "leave the site free, but it is kept aside because it is needed to "
         "place the grid and to validate the protocol. Waters, glycerol and "
         "sodium are crystallisation additives rather than parts of the enzyme."),
        ("Result",
         "Protein 9,184 atoms, FMN 124 atoms in four copies, XCJ 108 atoms in "
         "four copies, water 711 atoms. The receptor used for docking contains "
         "4,626 heavy atoms from chains A and B plus both FMN molecules. One "
         "residue, Lys254 of chain A, has an incomplete side chain in the "
         "crystal and was trimmed to alanine; it sits 12.6 &Aring; from the "
         "ligand and does not line the pocket."),
        ("For the write-up",
         "The crystal structure 7L00 was obtained from the Protein Data Bank. "
         "The protein and the essential cofactor FMN were retained, while the "
         "co-crystallised inhibitor XCJ, crystallographic water molecules, "
         "glycerol and sodium ions were removed during receptor preparation."),
    ]),
    (2, "Redocking validation", [
        ("Definitions",
         "<b>RMSD</b> = root mean square deviation, the average distance "
         "between corresponding atoms of two structures after they are "
         "compared in the same frame of reference. "
         "<b>Redocking</b> = removing the known inhibitor from a crystal "
         "structure and docking it back, then measuring how close the "
         "prediction comes to the experimental answer. "
         "<b>PDBQT</b> = the coordinate format AutoDock uses, which adds "
         "partial charges (Q) and atom types (T) to the PDB format."),
        ("What I did",
         "I generated a fresh three-dimensional conformer of XCJ before "
         "docking, so that it could not simply fall back into its crystal "
         "geometry, and docked it into the prepared receptor. I then measured "
         "the heavy-atom RMSD against the crystal pose without superposition, "
         "minimised over all symmetry-equivalent atom mappings so that a "
         "flipped phenyl ring is not counted as an error. Open Babel's "
         "<code>obrms</code> provided an independent second measurement."),
        ("Why",
         "A docking protocol that cannot reproduce a known answer cannot be "
         "trusted on unknown compounds. This is the single most important "
         "control in a docking study."),
        ("Result",
         "<b>RMSD = 0.82 &Aring;</b> for the top-ranked pose, with the second "
         "pose at 1.66 &Aring;. The conventional thresholds are: below 2 "
         "&Aring; excellent, 2 to 3 &Aring; acceptable, above 3 &Aring; "
         "re-optimise. The protocol passes comfortably. XCJ scored &minus;10.11 "
         "kcal/mol, which became the reference value for every later comparison."),
        ("For the write-up",
         "The co-crystallised inhibitor XCJ was redocked into the FabK active "
         "site to validate the docking protocol. The redocked pose reproduced "
         "the crystallographic binding conformation with an RMSD of 0.82 "
         "&Aring;, confirming the reliability of the docking method."),
    ]),
    (3, "Docking grid", [
        ("Definitions",
         "<b>Grid box</b> = the three-dimensional region the docking algorithm "
         "is allowed to search. A ligand cannot be placed outside it. "
         "<b>Grid centre</b> = the coordinates the box is built around."),
        ("What I did",
         "I centred the box on the geometric centre of the crystallographic "
         "XCJ, at (8.438, 15.490, 31.671). XCJ spans 11.9 &Aring; at its "
         "widest, so I first used a 22 &Aring; cube, adding 10 &Aring; of "
         "padding. When the compounds docked away from the cofactor I repeated "
         "the run with a 16 &Aring; cube, which forces the search into the "
         "FMN pocket itself. Both runs are reported."),
        ("Why",
         "Centring on the crystallographic ligand guarantees the search covers "
         "the real active site. Comparing a wide box with a narrow one "
         "distinguishes a compound that prefers a different region from one "
         "that cannot fit the pocket at all."),
        ("Result",
         "In the 22 &Aring; box the compounds settled about 9.7 &Aring; "
         "from the XCJ position and made no contact with FMN. In the 16 &Aring; "
         "box they occupied the true pocket, contacting 10 or 11 of the 22 "
         "residues that line the XCJ site, at a cost of 1.1 to 2.1 kcal/mol in "
         "docking score. The 16 &Aring; poses were used for all later analysis."),
        ("For the write-up",
         "Docking grids of 22 &Aring; and 16 &Aring; were centred on the "
         "position of the co-crystallised inhibitor XCJ to define the active "
         "site region for molecular docking studies."),
    ]),
    (4, "Molecular docking", [
        ("Definitions",
         "<b>Molecular docking</b> = a computational method that predicts the "
         "position, orientation and relative affinity of a small molecule "
         "within a binding site. <b>Docking score</b> = the estimated binding "
         "strength in kcal/mol, where a more negative value means tighter "
         "predicted binding. <b>Pose</b> = one predicted binding geometry. "
         "<b>Exhaustiveness</b> = how thoroughly the search algorithm samples. "
         "<b>MOL2</b> = Tripos molecular structure file, holding coordinates, "
         "atom types, bonds and charges; read by Discovery Studio, Chimera and "
         "MOE. <b>SDF</b> = structure data file. "
         "<b>TORSDOF</b> = torsional degrees of freedom, the number of "
         "rotatable bonds the docking treats as flexible."),
        ("What I did",
         "I prepared each compound from its MOL file: hydrogens added at pH "
         "7.4, a three-dimensional conformer generated, and the geometry "
         "minimised with the MMFF94 force field, then converted to PDBQT with "
         "Meeko. Docking used AutoDock Vina with exhaustiveness 32, 20 "
         "requested poses and a fixed random seed of 12345 so the run "
         "reproduces exactly. Every compound went through the same protocol "
         "that XCJ validated."),
        ("Why",
         "At pH 7.4 the benzoic acid of these compounds is deprotonated and "
         "carries a charge of &minus;1. Preparing them in that state matters, "
         "because the binding site turned out to be strongly anionic. A fixed "
         "seed makes the result reproducible rather than approximately "
         "repeatable."),
        ("Result",
         "Best scores in the pocket: 5sa23 &minus;7.99, 20sa23 &minus;7.30, "
         "17sa23 &minus;7.17 kcal/mol. The three differ by 0.8 kcal/mol, which "
         "is within Vina's normal uncertainty of 1 to 2 kcal/mol, so the "
         "docking score alone cannot separate them. This is why the study "
         "continued to MM/GBSA."),
        ("For the write-up",
         "Molecular docking was performed with AutoDock Vina using the "
         "validated protocol, with exhaustiveness 32 and twenty poses "
         "retained per compound. Binding poses were exported in MOL2 and PDB "
         "format for interaction analysis."),
    ]),
    (5, "MM/GBSA rescoring", [
        ("Definitions",
         "<b>MM/GBSA</b> = molecular mechanics with generalised Born surface "
         "area. <b>MM</b> = molecular mechanics, the force field description "
         "of the system. <b>GB</b> = generalised Born, an implicit model of "
         "water. <b>SA</b> = surface area, used for the non-polar part of "
         "solvation. <b>&Delta;G<sub>bind</sub></b> = G<sub>complex</sub> "
         "&minus; (G<sub>receptor</sub> + G<sub>ligand</sub>). "
         "<b>GAFF2</b> = general AMBER force field, version 2, used for small "
         "molecules. <b>ff14SB</b> = the AMBER protein force field. "
         "<b>AM1-BCC</b> = a semi-empirical method for assigning atomic "
         "partial charges. <b>igb=5</b> = the generalised Born variant used."),
        ("Energy terms",
         "<b>VDWAALS</b> = van der Waals, shape complementarity and buried "
         "surface. <b>EEL</b> = electrostatic interaction between charges. "
         "<b>EGB</b> = polar solvation, the cost of stripping water from "
         "charges on binding; it opposes EEL. <b>ESURF</b> = non-polar "
         "solvation, the burial of hydrophobic surface."),
        ("What I did",
         "I assigned GAFF2 parameters and AM1-BCC charges to each ligand and "
         "to FMN, built complex, receptor and ligand topologies with ff14SB, "
         "minimised each complex in implicit solvent with the protein backbone "
         "restrained, and ran single-structure MM/GBSA."),
        ("Why",
         "The Vina scoring function contains no electrostatic term at all. "
         "For anionic ligands in a charged pocket that is a serious omission, "
         "and MM/GBSA is the standard way to add it."),
        ("Result",
         "20sa23 &minus;32.02, 25sa23 &minus;28.53, 22sa23 &minus;23.35, "
         "5sa23 &minus;17.89, 17sa23 &minus;10.78 kcal/mol, against "
         "&minus;43.69 for the crystal inhibitor XCJ. The 21 kcal/mol spread "
         "separates the compounds where docking could not. "
         "Every EEL term for the compounds is positive, meaning "
         "electrostatics oppose binding, while XCJ's is favourable at "
         "&minus;42.51."),
        ("For the write-up",
         "MM/GBSA calculations were performed to estimate binding free "
         "energies and to refine the ranking obtained from molecular docking. "
         "Energies are relative ranking values: no entropy term is included, "
         "so they are not experimental binding free energies."),
    ]),
    (6, "Validation of the MM/GBSA result", [
        ("Why this was necessary",
         "The MM/GBSA ranking disagreed with the docking scores, and every "
         "electrostatic term came out positive. Before reporting a ranking "
         "that rests on an unexpected result, I tested three things that could "
         "have produced it artificially."),
        ("6A &mdash; control calculation",
         "I ran MM/GBSA on the redocked XCJ pose, which is known to be correct "
         "to 0.82 &Aring;. If a correct pose produced positive electrostatics "
         "too, the effect would belong to the setup rather than to the "
         "compounds. <b>Result:</b> XCJ gives EEL = &minus;42.51 kcal/mol, "
         "favourable. The positive values for the compounds are therefore "
         "real properties of those compounds."),
        ("6B &mdash; protonation state of His143",
         "<b>Protonation</b> = the addition of a hydrogen ion, H<sup>+</sup>. "
         "Histidine can carry the proton on either ring nitrogen or on both: "
         "<b>HID</b> is neutral with the proton on N&delta;, <b>HIE</b> is "
         "neutral with it on N&epsilon;, and <b>HIP</b> is doubly protonated "
         "and positively charged. His143 is FabK's catalytic residue, and the "
         "preparation had assumed the neutral form. Since the compounds are "
         "anions, a positive His143 could reverse the local electrostatics "
         "entirely. I therefore repeated every calculation in both states. "
         "<b>Result:</b> energies shift by 2 to 6 kcal/mol and the compound "
         "order is unchanged, so the ranking does not depend on this "
         "assumption."),
        ("6C &mdash; per-residue decomposition",
         "This divides the total binding energy among individual residues, "
         "showing which help and which oppose binding. I verified the output "
         "rather than accepting it: in this type of decomposition the "
         "contributions must sum back to the reported total. Van der Waals, "
         "electrostatics and polar solvation reconciled to within 0.014 "
         "kcal/mol across all eight calculations. The non-polar surface term "
         "does not decompose exactly, because with gbsa=2 the surface area is "
         "computed per atom by a route that is not additive per residue; that "
         "residual is bounded at 1.10 kcal/mol and reported separately."),
        ("Result",
         "The pocket is strongly anionic. Glu136 carries &minus;1 and the FMN "
         "phosphate carries &minus;2, and both repel the compounds' "
         "carboxylate: Glu136 contributes +12.7 to +34.1 kcal/mol of "
         "electrostatic repulsion and FMN +28.3 to +49.8, while for XCJ "
         "Glu136 is attractive at &minus;11.0. The one favourable anchor is "
         "His143, which reaches &minus;19.6 to &minus;26.3 kcal/mol when "
         "protonated, and it is not enough to overcome the rest."),
        ("For the write-up",
         "The binding free energy protocol was validated against the "
         "co-crystallised inhibitor, tested for sensitivity to the protonation "
         "state of the catalytic residue His143, and analysed by per-residue "
         "decomposition with an internal consistency check on the "
         "decomposition itself."),
    ]),
    (7, "Molecular dynamics", [
        ("Definitions",
         "<b>MD</b> = molecular dynamics, a simulation of how atoms move over "
         "time. Docking gives a single snapshot; MD gives a trajectory. "
         "<b>ns</b> = nanosecond, 10<sup>-9</sup> s; 50 ns is a reasonable "
         "length and 100 ns is a common publication standard. "
         "<b>RMSD</b> over a trajectory measures whether the complex stays "
         "stable. <b>RMSF</b> = root mean square fluctuation, the flexibility "
         "of each residue. Hydrogen bond occupancy measures how much of the "
         "simulation each bond persists for. Typical software: GROMACS, AMBER, "
         "NAMD, OpenMM."),
        ("Planned, not yet performed",
         "MD would test whether the predicted poses remain stable over time, "
         "and would allow binding free energies averaged over many "
         "conformations rather than one minimised structure. On the hardware "
         "used here it costs days of computation per compound. The next "
         "priority is extending the compound series, because the three-compound "
         "agreement reported below needs a larger set before it can be called "
         "a structure-activity relationship."),
    ]),
    (8, "Interaction diagrams", [
        ("Definitions",
         "<b>BIOVIA Discovery Studio Visualizer</b> = software that produces "
         "two-dimensional interaction maps of a protein-ligand complex, "
         "showing hydrogen bonds, hydrophobic contacts, &pi;-&pi; stacking, "
         "&pi;-cation interactions and van der Waals contacts."),
        ("Files prepared",
         "The receptor <code>7L00_AB_protein_FMN.pdb</code> opens together "
         "with the best pose of each of the five compounds, "
         "<code>docking_focused_&lt;compound&gt;/&lt;compound&gt;_best.mol2</code>. "
         "The minimised complexes used for the MM/GBSA energies are also "
         "available as single PDB files, "
         "<code>mmgbsa_focused/&lt;compound&gt;/complex_min.pdb</code>, each "
         "containing protein, FMN and ligand together."),
        ("Planned, not yet performed",
         "The input files are prepared and listed above; the diagrams "
         "themselves have not yet been generated. They will provide the visual "
         "counterpart to the per-residue energies in section 6, which place "
         "the carboxylate of 5sa23 and 17sa23 closest to Glu136 and the FMN "
         "phosphate, and that of 20sa23 and 25sa23 furthest from both."),
        ("Note",
         "No sentence is offered here for the manuscript, because the "
         "diagrams do not yet exist. One will be written when they do."),
    ]),
    (9, "Comparison with measured antibacterial activity", [
        ("Definitions",
         "<b>MIC</b> = minimum inhibitory concentration, the lowest "
         "concentration of a compound that prevents visible bacterial growth, "
         "reported in &micro;g/mL. A lower MIC means a more potent compound. "
         "As a rough guide: below 1 &micro;g/mL excellent, 1 to 10 good, 10 to "
         "50 moderate, above 50 weak. "
         "<b>Pearson r</b> measures linear correlation; <b>Spearman r</b> "
         "measures whether the rank order agrees. Correlation is computed "
         "against log<sub>10</sub>(MIC) because concentration spans orders of "
         "magnitude."),
        ("What I did",
         "I compared every prediction the study produced against the MIC "
         "values reported in the source paper for these five compounds "
         "against <i>Staphylococcus aureus</i> ATCC 29213."),
        ("Result",
         "The MM/GBSA ranking reproduces the experimental order exactly. The "
         "docking score does not, and the count of close polar contacts gives "
         "the order backwards. Full comparison in the next section."),
        ("For the write-up",
         "The computational results were compared with experimentally "
         "determined MIC values to evaluate whether predicted binding "
         "affinities corresponded to antibacterial potency."),
    ]),
]

GLOSSARY = [
    ("&Aring;", "&aring;ngstr&ouml;m, 10<sup>-10</sup> m"),
    ("AM1-BCC", "semi-empirical atomic partial charge method"),
    ("EEL", "electrostatic energy term"),
    ("EGB", "polar solvation energy term (generalised Born)"),
    ("ESURF", "non-polar solvation energy term (surface area)"),
    ("FabK", "enoyl-ACP reductase II, a bacterial fatty acid biosynthesis enzyme"),
    ("ff14SB", "AMBER protein force field"),
    ("FMN", "flavin mononucleotide, the FabK cofactor"),
    ("GAFF2", "general AMBER force field for small molecules"),
    ("HID / HIE / HIP", "histidine protonated on N&delta; / on N&epsilon; / on both (charged)"),
    ("kcal/mol", "kilocalorie per mole, the energy unit used throughout"),
    ("MD", "molecular dynamics"),
    ("MIC", "minimum inhibitory concentration, &micro;g/mL"),
    ("MM/GBSA", "molecular mechanics / generalised Born surface area"),
    ("MOL2", "Tripos molecular structure file"),
    ("ns", "nanosecond, 10<sup>-9</sup> s"),
    ("PDB", "Protein Data Bank, and its coordinate file format"),
    ("PDBQT", "AutoDock coordinate format with charges and atom types"),
    ("RMSD", "root mean square deviation"),
    ("RMSF", "root mean square fluctuation"),
    ("TORSDOF", "torsional degrees of freedom (rotatable bonds)"),
    ("VDWAALS", "van der Waals energy term"),
    ("XCJ", "the co-crystallised FabK inhibitor in 7L00"),
]

# ============================================================
# Files to gather
# ============================================================
SCRIPTS = [
    ("prepare_7L00_structure.sh", "1", "Downloads 7L00 and separates protein, FMN, XCJ, water and other components."),
    ("step2_7L00_receptor_redock.sh", "2-3", "Builds the receptor, characterises the site, sets the grid, redocks XCJ and measures RMSD."),
    ("step3_7L00_dock_compounds.sh", "4", "Prepares and docks the compounds, exports poses and complexes, analyses contacts, ranks and graphs."),
    ("step3b_pose_scan.py", "4", "Checks anchor-residue contacts across all twenty poses rather than the best one alone."),
    ("step4_mmgbsa.sh", "5", "Parameterises ligands and FMN, builds topologies, minimises and runs MM/GBSA."),
    ("step6_mmgbsa_checks.sh", "6", "Control calculation on XCJ, both His143 protonation states, and decomposition."),
    ("step6c_decomposition.sh", "6", "Per-residue decomposition with a consistency check on the decomposition itself."),
    ("diagnose_decomp.sh", "6", "Read-only diagnostic used to locate the cause of a decomposition failure."),
    ("diagnose_decomp_residual.py", "6", "Read-only component-wise reconciliation of decomposition sums against totals."),
    ("item9_mic_correlation.py", "9", "Compares every prediction against the measured MIC values."),
    ("build_report.sh", "-", "Builds this page."),
]

# Superseded, redundant or scratch output that is not carried into the report.
SKIP_PREFIXES = ("00_",)
SKIP_EXACT = {"decomp_diagnostic.txt"}

SECTIONS = [
    ("Structure preparation", ["01_", "02_", "03_", "04_", "05_", "06_", "07_"]),
    ("Receptor and binding site", ["08_", "09_", "10_", "11_", "12_", "14_"]),
    ("Redocking validation", ["13_", "15_", "16_", "17_", "17b_", "18_"]),
    ("Compound preparation", ["19_", "20_", "21_"]),
    ("Docking runs", ["22_", "23_", "24_"]),
    ("Contacts and ranking", ["25_", "26_", "27_", "34_"]),
    ("MM/GBSA and validation", ["35_", "36_", "37_", "38_", "39_", "40_", "41_"]),
    ("MIC comparison", ["42_"]),
]

root = Path(".")
txts = sorted(
    p for p in root.rglob("*.txt")
    if p.is_file() and "txt_outputs" in p.parts and "focused" in p.parts
    and p.stat().st_size > 0
    and not p.name.startswith(SKIP_PREFIXES) and p.name not in SKIP_EXACT
)
pngs = sorted(p for p in root.rglob("*.png")
              if p.is_file() and p.parts[0] == "graphs_focused")


def run_label(path: Path) -> str:
    return ""


def section_of(path: Path) -> str:
    for title, prefixes in SECTIONS:
        if any(path.name.startswith(p) for p in prefixes):
            return title
    return "Other output"


def strip_identity(text: str) -> str:
    lines = text.splitlines()
    rule = lambda l: set(l.strip()) == {"="} and len(l.strip()) > 20
    i = 0
    while i < len(lines):
        if rule(lines[i]):
            for j in range(i + 1, min(i + 9, len(lines))):
                if rule(lines[j]):
                    if any("siam" in lines[k].lower() for k in range(i + 1, j)):
                        del lines[i:j + 1]
                        i -= 1
                    break
        i += 1
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return "\n".join(lines)


grouped = {}
for p in txts:
    grouped.setdefault(section_of(p), []).append(p)

png_groups = {}
for p in pngs:
    png_groups.setdefault(p.parts[0], []).append(p)

# ============================================================
# HTML
# ============================================================
E = html.escape
P = []

P.append(f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Molecular docking study: 7L00 FabK &mdash; {E(AUTHOR)}</title>
<style>
  :root {{ --ink:#1a1a1a; --muted:#666; --line:#dcdcdc; --bg:#fff;
           --card:#fafafa; --accent:#1f6f4a; --warn:#b4432b; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--ink);
          font:15.5px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
  .wrap {{ max-width:1080px; margin:0 auto; padding:0 26px 90px; }}
  header {{ border-bottom:3px solid var(--ink); padding:40px 0 22px; }}
  h1 {{ margin:0 0 4px; font-size:29px; letter-spacing:-.015em; line-height:1.25; }}
  .sub {{ font-size:16.5px; color:#333; margin:14px 0 20px; padding:13px 17px;
           background:#eef6f1; border-left:5px solid var(--accent);
           border-radius:0 7px 7px 0; }}
  .sub b {{ color:var(--accent); }}
  .who {{ font-size:14.5px; line-height:1.55; }}
  .who .name {{ font-size:16px; font-weight:700; }}
  .who .aff {{ color:var(--muted); }}
  .who .aff a {{ color:var(--accent); text-decoration:none; }}
  .who .aff a:hover {{ text-decoration:underline; }}
  .who .links {{ margin:5px 0 7px; font-size:13.5px; }}
  .who .links a {{ color:var(--accent); text-decoration:none; font-weight:600;
                   margin-right:16px; }}
  .who .links a:hover {{ text-decoration:underline; }}
  footer a {{ color:var(--accent); }}
  .meta {{ margin-top:14px; font-size:13.5px; color:var(--muted); }}
  .meta a {{ color:var(--accent); }}
  nav {{ position:sticky; top:0; background:rgba(255,255,255,.97);
         backdrop-filter:blur(6px); padding:12px 0; border-bottom:1px solid var(--line);
         margin-bottom:10px; z-index:20; }}
  nav a {{ display:inline-block; margin:3px 15px 3px 0; color:var(--accent);
           text-decoration:none; font-size:13px; font-weight:600; }}
  nav a:hover {{ text-decoration:underline; }}
  h2 {{ margin:52px 0 6px; font-size:22px; padding-bottom:9px;
        border-bottom:2px solid var(--line); }}
  h3 {{ margin:30px 0 8px; font-size:17px; }}
  .step {{ border:1px solid var(--line); border-radius:8px; padding:20px 22px;
           margin:18px 0; background:var(--card); }}
  .step h3 {{ margin:0 0 14px; font-size:18.5px; display:flex; gap:11px;
              align-items:baseline; }}
  .num {{ background:var(--ink); color:#fff; border-radius:5px; font-size:13px;
          padding:2px 9px; font-weight:700; flex:none; }}
  .pending {{ background:#f3ead6; color:#7a5a12; border:1px solid #e0cfa4;
              border-radius:4px; font-size:11.5px; font-weight:700; padding:2px 9px;
              text-transform:uppercase; letter-spacing:.04em; margin-left:auto;
              flex:none; }}
  .pending-step {{ background:#fdfbf5; border-style:dashed; }}
  .block {{ margin:13px 0; }}
  .block .lbl {{ font-weight:700; font-size:13px; text-transform:uppercase;
                 letter-spacing:.045em; color:var(--muted); margin-bottom:3px; }}
  table {{ border-collapse:collapse; width:100%; margin:16px 0; font-size:14px; }}
  th, td {{ border:1px solid var(--line); padding:8px 11px; text-align:left; }}
  th {{ background:#f0f0f0; font-size:13px; }}
  td.n {{ text-align:right; font-variant-numeric:tabular-nums; }}
  tr.best td {{ background:#eaf5ee; }}
  tr.worst td {{ background:#fbeeea; }}
  tr.ref td {{ background:#f4f4f4; color:#555; font-style:italic; }}
  .key {{ border-left:5px solid var(--accent); background:#eef6f1;
          padding:18px 22px; border-radius:0 8px 8px 0; margin:24px 0; }}
  .key h3 {{ margin:0 0 9px; color:var(--accent); font-size:18px; }}
  .alert {{ border:2px solid var(--warn); background:#fdf2ef; padding:18px 22px;
             border-radius:8px; margin:26px 0 8px; }}
  .alert h3 {{ margin:0 0 10px; color:var(--warn); font-size:17.5px; }}
  .alert p {{ margin:9px 0; font-size:14.5px; }}
  .caveat {{ border-left:5px solid var(--warn); background:#fdf2ef;
             padding:15px 20px; border-radius:0 8px 8px 0; margin:22px 0;
             font-size:14.5px; }}
  details {{ border:1px solid var(--line); border-radius:6px; margin:9px 0;
             background:var(--card); }}
  details[open] {{ background:#fff; }}
  summary {{ cursor:pointer; padding:11px 14px; font-weight:600; font-size:14px;
             list-style:none; display:flex; justify-content:space-between; gap:14px; }}
  summary::-webkit-details-marker {{ display:none; }}
  summary::before {{ content:"\\25B8"; margin-right:9px; color:var(--muted); }}
  details[open] summary::before {{ content:"\\25BE"; }}
  .path {{ color:var(--muted); font-weight:400; font-size:11.5px;
           font-family:ui-monospace,Menlo,Consolas,monospace; }}
  pre {{ margin:0; padding:15px 17px; overflow-x:auto; font-size:12.5px;
         line-height:1.5; font-family:ui-monospace,Menlo,Consolas,monospace;
         border-top:1px solid var(--line); background:#fff; white-space:pre; }}
  code {{ font-family:ui-monospace,Menlo,Consolas,monospace; font-size:13px;
          background:#f0f0f0; padding:1px 5px; border-radius:3px; }}
  figure {{ margin:0 0 26px; }}
  figure img {{ width:100%; border:1px solid var(--line); border-radius:6px; }}
  figcaption {{ color:var(--muted); font-size:13px; margin-top:7px; }}
  .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(430px,1fr));
           gap:26px; }}
  .tag {{ display:inline-block; background:#e8eef5; color:#2c4a68;
          border-radius:4px; padding:1px 7px; font-size:11px;
          font-weight:600; margin-left:8px; }}
  .bar {{ display:flex; gap:10px; margin:20px 0 4px; flex-wrap:wrap; }}
  button {{ font:inherit; font-size:13px; padding:7px 15px; border:1px solid var(--line);
            background:#fff; border-radius:6px; cursor:pointer; }}
  button:hover {{ background:var(--card); }}
  dl.gloss {{ display:grid; grid-template-columns:max-content 1fr; gap:6px 20px;
              font-size:14px; margin:14px 0; }}
  dl.gloss dt {{ font-weight:700; font-family:ui-monospace,Menlo,Consolas,monospace; }}
  dl.gloss dd {{ margin:0; color:#333; }}
  footer {{ margin-top:60px; padding-top:20px; border-top:2px solid var(--line);
            color:var(--muted); font-size:13px; }}
  @media print {{
    nav, .bar {{ display:none; }}
    .step, details, table, figure {{ break-inside:avoid; }}
    details:not([open]) > *:not(summary) {{ display:block; }}
  }}
</style></head><body><div class="wrap">

<header>
<h1>Molecular docking and binding free energy study of pyrazole benzoic acid
derivatives against <i>Clostridioides difficile</i> FabK (PDB 7L00)</h1>
<div class="sub">Computational study of fatty acid biosynthesis inhibition at
FabK (PDB 7L00). Across five compounds, MM/GBSA reproduced the measured
antibacterial potency order exactly, with no pair transposed
(Spearman&nbsp;1.00, Pearson&nbsp;0.99):<br>
measured <b>20sa23 (MIC 1&nbsp;&micro;g/mL) &gt; 25sa23 (MIC
2&nbsp;&micro;g/mL) &gt; 22sa23 (MIC 4&nbsp;&micro;g/mL) &gt; 5sa23 (MIC
16&nbsp;&micro;g/mL) &gt; 17sa23 (MIC
32&nbsp;&micro;g/mL)</b>, predicted
<b>20sa23 &gt; 25sa23 &gt; 22sa23 &gt; 5sa23 &gt; 17sa23</b>.
The docking score alone gave no relationship (Pearson&nbsp;0.20).</div>
<div class="who">
  <div class="name">{E(AUTHOR)}</div>
  <div class="links">
    <a href="https://{AUTHOR_SITE}">{AUTHOR_SITE}</a>
    <a href="{AUTHOR_GITHUB}">GitHub</a>
    <a href="{AUTHOR_ORCID}">ORCID 0000-0001-5216-6879</a>
  </div>
  <div class="aff">{'<br>'.join(AFFIL)}</div>
</div>
<div class="meta">
  Report compiled {E(NOW)}<br>
  Reference paper: {E(PAPER_TITLE)}<br>
  <a href="{PAPER_DOI}">{PAPER_DOI}</a>
</div>
</header>

<nav>
  <a href="#summary">Summary</a>
  <a href="#methods">Methods</a>
  <a href="#findings">Key findings</a>
  <a href="#graphs">Graphs</a>
  <a href="#scripts">Scripts</a>
  <a href="#output">Output files</a>
  <a href="#glossary">Glossary</a>
</nav>

<div class="bar">
  <button onclick="document.querySelectorAll('details').forEach(d=>d.open=true)">Expand all</button>
  <button onclick="document.querySelectorAll('details').forEach(d=>d.open=false)">Collapse all</button>
  <button onclick="window.print()">Print / save as PDF</button>
</div>
""")

# ---------------- summary ----------------
P.append("""
<div class="alert">
<h3>Target and organism caveat &mdash; read first</h3>
<p>FabK belongs to a family of closely related enoyl-ACP reductases found
across many bacterial species. The structure used here is
<i>Clostridioides difficile</i> FabK (PDB 7L00), while the MIC values the
predictions are compared against were measured in <i>Staphylococcus aureus</i>
ATCC&nbsp;29213. A compound may bind one homologue differently from another,
and whole-cell potency in one organism need not reflect binding to an enzyme
from a different one.</p>
<p>The comparison is therefore between a binding model for one protein and
antibacterial activity in a different organism. Properly matched work would
use the enzyme of the organism the MIC was measured in, or MIC data for the
organism the structure came from. The agreement reported below should be read
with that mismatch in mind.</p>
</div>

<h2 id="summary">Summary</h2>
<p>Three pyrazole benzoic acid derivatives (5sa23, 17sa23 and 20sa23) were
docked against FabK, the enoyl-ACP reductase of bacterial fatty acid
biosynthesis, using the crystal structure 7L00. The docking protocol was
validated by redocking the co-crystallised inhibitor XCJ, which reproduced the
experimental pose to 0.82&nbsp;&Aring; RMSD. The compounds were then rescored
with MM/GBSA, and that result was itself validated with a control calculation,
a test of the catalytic histidine's protonation state, and a per-residue energy
decomposition with an internal consistency check.</p>

<p><b>Scope.</b> Sections 1 to 6 and section 9 were carried out. Molecular
dynamics (section 7) and the Discovery Studio interaction diagrams (section 8)
are planned and have not yet been performed; they are marked as such below.</p>

<p>The docking scores could not separate the compounds: they fell within
1.2&nbsp;kcal/mol of one another, inside the method's own uncertainty.
MM/GBSA separated them by 21&nbsp;kcal/mol and produced an order that matches
the measured antibacterial potency exactly. The reason is electrostatic. The
FabK pocket is strongly anionic, carrying Glu136 at &minus;1 and the FMN
phosphate at &minus;2, while these compounds are benzoate anions. The most
potent compound is the one that keeps its carboxylate furthest from that
negative charge.</p>
""")

# ---------------- methods ----------------
P.append('<h2 id="methods">Methods</h2>')
for num, title, blocks in METHODS:
    flag = ('<span class="pending">not yet performed</span>'
            if num in PENDING else "")
    cls = "step pending-step" if num in PENDING else "step"
    P.append(f'<div class="{cls}"><h3><span class="num">{num}</span>'
             f'<span>{title}</span>{flag}</h3>')
    for lbl, body in blocks:
        P.append(f'<div class="block"><div class="lbl">{lbl}</div>{body}</div>')
    P.append("</div>")

# ---------------- key findings ----------------
rows = []
for i, (name, sub, mic, vina, gbsa, eel, glu, fmn, polar) in enumerate(COMPOUNDS):
    cls = "best" if i == 0 else ("worst" if i == len(COMPOUNDS) - 1 else "")
    verdict = ("<b>Most potent</b>" if i == 0 else
               ("<b>Least potent</b>" if i == len(COMPOUNDS) - 1 else ""))
    rows.append(
        f'<tr class="{cls}"><td><b>{name}</b></td><td>{sub}</td>'
        f'<td class="n">{mic}</td>'
        f'<td class="n">{vina:.2f}</td><td class="n">{gbsa:.2f}</td>'
        f'<td class="n">{eel:+.1f}</td><td class="n">{glu:+.1f}</td>'
        f'<td class="n">{fmn:+.1f}</td><td class="n">{polar}</td>'
        f'<td>{verdict}</td></tr>')
n, sub, mic, vina, gbsa, eel, glu, fmn, polar = XCJ
rows.append(
    f'<tr class="ref"><td>{n}</td><td>{sub}</td><td class="n">{mic}</td>'
    f'<td class="n">{vina:.2f}</td><td class="n">{gbsa:.2f}</td>'
    f'<td class="n">{eel:+.1f}</td><td class="n">{glu:+.1f}</td>'
    f'<td class="n">{fmn:+.1f}</td><td class="n">{polar}</td>'
    f'<td>Validation reference</td></tr>')

P.append(f"""
<h2 id="findings">Key findings</h2>

<div class="key">
<h3>MM/GBSA reproduces the experimental potency order; the docking score does not</h3>
<p>Measured against <i>S. aureus</i> ATCC&nbsp;29213, as reported in the source
paper:<br>
<b>20sa23 (MIC 1&nbsp;&micro;g/mL) &gt; 25sa23 (MIC 2&nbsp;&micro;g/mL) &gt;
22sa23 (MIC 4&nbsp;&micro;g/mL) &gt; 5sa23
(MIC 16&nbsp;&micro;g/mL) &gt; 17sa23 (MIC 32&nbsp;&micro;g/mL)</b></p>
<p>Predicted by MM/GBSA from the structures alone:<br>
<b>20sa23 (&minus;32.0) &gt; 25sa23 (&minus;28.5) &gt; 22sa23 (&minus;23.4)
&gt; 5sa23 (&minus;17.9) &gt; 17sa23
(&minus;10.8&nbsp;kcal/mol)</b></p>
<p>All five positions are correct, from the most potent compound to the least
potent, across the full 32-fold range of measured MIC. No adjacent pair is
transposed. Spearman 1.00, Pearson 0.99 against
log<sub>10</sub>(MIC).</p>
<p>The docking score gave no relationship (Pearson 0.20), and the count of
close polar contacts pointed the wrong way (&minus;0.87).</p>
</div>

<div class="key">
<h3>The isomer pair: same atoms, different potency, predicted correctly</h3>
<p>22sa23 and 25sa23 share the molecular formula
C<sub>24</sub>H<sub>16</sub>F<sub>5</sub>N<sub>3</sub>O<sub>2</sub> and differ
only in whether the aniline fluorine sits at position 4 or position 5 relative
to the CF<sub>3</sub> group. Moving that one atom changes the Glu136
electrostatic repulsion from <b>+32.5 to +14.2&nbsp;kcal/mol</b> and the
binding energy from &minus;23.4 to &minus;28.5. The measured MICs are 4 and
2&nbsp;&micro;g/mL respectively, in the predicted direction.</p>
<p>Because the two molecules are otherwise identical &mdash; same atoms, same
charge, same rotatable bonds &mdash; the difference can only be geometric, which
makes this the cleanest single test in the study.</p>
</div>

<h3>Final comparison</h3>
<table>
<tr><th>Compound</th><th>Aniline</th><th>MIC<br>&micro;g/mL</th>
<th>Vina<br>kcal/mol</th>
<th>MM/GBSA<br>kcal/mol</th><th>EEL<br>kcal/mol</th>
<th>Glu136<br>EEL</th><th>FMN<br>EEL</th><th>Polar<br>contacts</th>
<th></th></tr>
{''.join(rows)}
</table>
<p style="font-size:14px;color:#555">Sorted by measured potency. Lower MIC means
more potent; more negative energy means tighter predicted binding; positive EEL
means electrostatic repulsion. Per-residue values are from the neutral His143
calculation. All five compounds share the same core, a 4-fluorophenyl pyrazole
linked to benzoic acid, and differ only in the aniline substitution. Both
3,5-disubstituted compounds sit near +13&nbsp;kcal/mol of Glu136 repulsion;
every other compound is above +32.</p>

<h3>Which prediction tracked potency</h3>
<table>
<tr><th>Prediction</th><th>Order produced</th>
<th>Pearson r<br>vs log(MIC)</th><th>Spearman</th></tr>
<tr><td><b>MM/GBSA total energy</b></td>
<td>20 &gt; 25 &gt; 22 &gt; 5 &gt; 17</td>
<td class="n"><b>0.99</b></td><td class="n"><b>1.00</b></td></tr>
<tr><td>FMN electrostatic contribution</td>
<td>20 &gt; 25 &gt; 22 &gt; 5 &gt; 17</td>
<td class="n">0.99</td><td class="n">1.00</td></tr>
<tr><td>MM/GBSA electrostatic term</td>
<td>20 &gt; 25 &gt; 22 &gt; 5 &gt; 17</td>
<td class="n">0.93</td><td class="n">1.00</td></tr>
<tr><td>His143 electrostatic contribution</td>
<td>20 &gt; 25 &gt; 22 &gt; 5 &gt; 17</td>
<td class="n">0.91</td><td class="n">1.00</td></tr>
<tr><td>Glu136 electrostatic contribution</td>
<td>20 &gt; 25 &gt; 22 &gt; 17 &gt; 5</td>
<td class="n">0.85</td><td class="n">0.90</td></tr>
<tr><td>Vina docking score</td>
<td>22 &gt; 25 &gt; 5 &gt; 20 &gt; 17</td>
<td class="n">0.20</td><td class="n">0.30</td></tr>
<tr><td>Van der Waals term</td>
<td>22 &gt; 5 &gt; 17 &gt; 25 &gt; 20</td>
<td class="n">&minus;0.44</td><td class="n">&minus;0.60</td></tr>
<tr class="worst"><td>Count of close polar contacts</td>
<td>17 &gt; 22 &gt; 5 &gt; 20 &gt; 25</td>
<td class="n">&minus;0.87</td><td class="n">&minus;0.87</td></tr>
</table>
<p style="font-size:14px;color:#555">A working predictor gives a positive r:
lower MIC together with more negative energy. Compound numbers abbreviated
(20 = 20sa23, and so on).</p>

<h3>MM/GBSA energy components</h3>
<table>
<tr><th>Compound</th><th>van der Waals<br>VDWAALS</th>
<th>Electrostatic<br>EEL</th><th>Polar solvation<br>EGB</th>
<th>Non-polar<br>ESURF</th><th>Total<br>&Delta;G<sub>bind</sub></th></tr>
{''.join(
    f'<tr class="{"ref" if n == "XCJ" else ""}">'
    f'<td>{"<b>" + n + "</b>" if n != "XCJ" else n}</td>'
    f'<td class="n">{v:.2f}</td><td class="n">{e:+.2f}</td>'
    f'<td class="n">{g:+.2f}</td><td class="n">{su:.2f}</td>'
    f'<td class="n"><b>{t:.2f}</b></td></tr>'
    for n, v, e, g, su, t in COMPONENTS)}
</table>
<p style="font-size:14px;color:#555">All values kcal/mol, neutral His143, sorted
by measured potency. EEL is positive for every compound, meaning electrostatics
oppose binding, and it is the term that separates them: it spans 55 to 127
across the series while van der Waals stays within 14. For the crystal
inhibitor XCJ, EEL is favourable at &minus;42.51, which is what confirms the
positive values are a property of these compounds rather than of the
calculation. EGB opposes EEL because charges must shed their water shell to
bind, so the two partly cancel and the total is the difference between two
large numbers.</p>

<h3>Per-residue electrostatic contributions</h3>
<table>
<tr><th>Compound</th><th>Glu136</th><th>His143<br>neutral</th>
<th>His143<br>protonated</th><th>FMN</th></tr>
{''.join(
    f'<tr class="{"ref" if n == "XCJ" else ""}">'
    f'<td>{"<b>" + n + "</b>" if n != "XCJ" else n}</td>'
    f'<td class="n">{glu:+.2f}</td><td class="n">{h_e:+.2f}</td>'
    f'<td class="n">{h_p:+.2f}</td><td class="n">{fmn:+.2f}</td></tr>'
    for n, glu, h_e, h_p, fmn in RESIDUES)}
</table>
<p style="font-size:14px;color:#555">Electrostatic contribution of each residue
to the binding energy, kcal/mol. Positive opposes binding. Glu136 carries
&minus;1 and the FMN phosphate &minus;2, and both repel the compounds'
carboxylate; for XCJ, which has no carboxylate in that position, Glu136 is
attractive instead. His143 is the one favourable anchor, and only when
protonated and positively charged, where it reaches &minus;20 to &minus;29 for
the compounds against &minus;7 for XCJ. These values come from the
decomposition that was verified to reconcile with the totals to within 0.014
kcal/mol.</p>

<h3>Sensitivity to the His143 protonation state</h3>
<table>
<tr><th>Compound</th><th>His143 neutral<br>(HIE)</th>
<th>His143 protonated<br>(HIP)</th><th>Difference</th></tr>
{''.join(
    f'<tr class="{"ref" if n == "XCJ" else ""}">'
    f'<td>{"<b>" + n + "</b>" if n != "XCJ" else n}</td>'
    f'<td class="n">{hie:.2f}</td><td class="n">{hip:.2f}</td>'
    f'<td class="n">{hip - hie:+.2f}</td></tr>'
    for n, hie, hip in PROTONATION)}
</table>
<p style="font-size:14px;color:#555">Total binding energy, kcal/mol, with the
catalytic histidine modelled as neutral and as protonated. Every compound binds
more tightly with His143 protonated, by 2 to 6 kcal/mol, which is expected
since a positively charged histidine attracts an anion. The ordering of the
compounds is unchanged between the two states, so the ranking does not depend
on this assumption. The preparation assigns the neutral form by default, and
testing the alternative was one of the three checks applied before the MM/GBSA
result was accepted.</p>

<h3>Why the contact count points the wrong way</h3>
<p>Counting contacts by distance alone treats any close nitrogen-oxygen pair as
a favourable hydrogen bond. At pH&nbsp;7.4 the compounds' benzoate carries
&minus;1 and Glu136 carries &minus;1, so a close approach between them is
electrostatic repulsion, not a hydrogen bond. Across all five compounds the
relationship is inverted (r = &minus;0.87): 17sa23 makes the most such contacts
and is the least potent, while 20sa23 and 25sa23 make the fewest and are the
two most potent. The AutoDock Vina scoring function contains no electrostatic term
either, which is why it gave no useful signal. MM/GBSA models charge
explicitly, and that is the difference.</p>

<div class="caveat">
<b>Target and organism.</b> The caveat at the top of this page applies to every
number in this section: the structure is <i>C. difficile</i> FabK and the MIC
values are from <i>S. aureus</i>. The correlation is between a binding model
for one enzyme and whole-cell activity in an organism that may not depend on
that enzyme in the same way.
</div>

<div class="caveat">
<b>How far this can be taken.</b> Five compounds is still a small set, and a
correlation at this size carries limited statistical weight. What supports the
result is that the ordering follows from a specific mechanism &mdash; an anionic
pocket rejecting an anionic ligand, with the aniline substitution pattern
determining how close the carboxylate comes to Glu136 and the FMN phosphate
&mdash; and that this mechanism predicted the isomer pair correctly. The
compounds were analysed and ranked before their measured values were
consulted. It is also worth noting that all five compounds remain
electrostatically repelled by
this pocket relative to XCJ, so whether they inhibit FabK by binding here
remains open. The source paper assigned the fatty acid biosynthesis pathway by
CRISPRi, which identifies the pathway rather than the individual enzyme.
Finally, MIC reflects whole-cell activity, which depends on membrane
permeability, efflux and metabolic stability as well as target binding.
</div>
""")

# ---------------- graphs ----------------
if pngs:
    P.append('<h2 id="graphs">Graphs</h2>')
    for folder in sorted(png_groups):
        P.append(f'<h3>{E(folder)}/</h3>')
        P.append('<div class="grid">')
        for p in png_groups[folder]:
            b64 = base64.b64encode(p.read_bytes()).decode()
            P.append(f'<figure><img src="data:image/png;base64,{b64}" '
                     f'alt="{E(p.name)}"><figcaption>{E(p.name)}</figcaption></figure>')
        P.append("</div>")

# ---------------- scripts ----------------
P.append('<h2 id="scripts">Scripts</h2>'
         '<p>Every stage of this study was run from scripts rather than a '
         'graphical interface, so the whole pipeline is reproducible from the '
         'source structure. The full source of each is below.</p>')
found_any = False
for fname, step, desc in SCRIPTS:
    f = Path(fname)
    if not f.is_file():
        continue
    found_any = True
    step_tag = f'<span class="tag">step {step}</span>' if step != "-" else ""
    P.append(
        f'<details><summary><span>{E(fname)}{step_tag}</span>'
        f'<span class="path">{len(f.read_text(errors="replace").splitlines())} lines</span></summary>'
        f'<div style="padding:12px 16px 0;font-size:14px;color:#555">{desc}</div>'
        f'<pre>{E(f.read_text(errors="replace"))}</pre></details>')
if not found_any:
    P.append('<p style="color:#888">No script files found in this folder.</p>')

# ---------------- text output ----------------
P.append('<h2 id="output">Output files</h2>'
         '<p>Each stage writes a plain-text report recording what it did, what '
         'it produced and any warnings, so that every number above can be '
         'traced back to the run that produced it.</p>')
for title, _ in SECTIONS + [("Other output", [])]:
    if title not in grouped:
        continue
    P.append(f'<h3>{E(title)}</h3>')
    for p in sorted(grouped[title]):
        body = strip_identity(p.read_text(errors="replace"))
        lbl = run_label(p)
        tag = f'<span class="tag">{lbl}</span>' if lbl else ""
        P.append(f'<details><summary><span>{E(p.name)}{tag}</span>'
                 f'<span class="path">{E(str(p))}</span></summary>'
                 f'<pre>{E(body)}</pre></details>')

# ---------------- glossary ----------------
P.append('<h2 id="glossary">Glossary</h2><dl class="gloss">')
for term, meaning in GLOSSARY:
    P.append(f'<dt>{term}</dt><dd>{meaning}</dd>')
P.append("</dl>")

P.append(f"""
<footer>
<b>{E(AUTHOR)}</b> &middot;
<a href="https://{AUTHOR_SITE}">{AUTHOR_SITE}</a> &middot;
<a href="{AUTHOR_GITHUB}">GitHub</a> &middot;
<a href="{AUTHOR_ORCID}">ORCID 0000-0001-5216-6879</a><br>
{'<br>'.join(AFFIL)}<br><br>
{len(txts)} text reports and {len(pngs)} graphs, collected from their original
locations; no file was moved or altered.<br>
Docking and MM/GBSA values are computational scores, not experimentally
measured binding free energies.
</footer>
</div></body></html>""")

out_path.write_text("".join(P), encoding="utf-8")
size = out_path.stat().st_size / 1024
print(f"{len(txts)} reports, {len(pngs)} graphs, "
      f"{sum(1 for f, _, _ in SCRIPTS if Path(f).is_file())} scripts "
      f"-> {out_path} ({size:.0f} KB)")
PY

echo
echo "Report ready: $(pwd)/$OUT"

if (( OPEN )); then
    if command -v wslview >/dev/null 2>&1; then
        wslview "$OUT" 2>/dev/null &
    elif command -v explorer.exe >/dev/null 2>&1; then
        explorer.exe "$(wslpath -w "$OUT")" 2>/dev/null &
    elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open "$OUT" 2>/dev/null &
    fi
    sleep 1
fi
