# FabK docking study (PDB 7L00)

Molecular docking and MM/GBSA study of five pyrazole benzoic acid derivatives
against *Clostridioides difficile* FabK, an enzyme of bacterial fatty acid
biosynthesis.

**Full report:** https://siamchowdhury-1.github.io/fabk-docking-7L00/

The docking protocol was validated by redocking the co-crystallised inhibitor
XCJ, reproducing the experimental pose to 0.82 Å RMSD. MM/GBSA rescoring
reproduced the measured antibacterial potency order of all five compounds
exactly, across a 32-fold range of MIC; the docking score alone did not.

The pocket is strongly anionic, carrying Glu136 at −1 and the FMN phosphate
at −2, while these compounds are benzoate anions. Every compound is
electrostatically repelled, and how strongly depends on where the aniline
substitution places the carboxylate. That is what separates them, and it is
why a scoring function without an explicit electrostatic term could not.

| Compound | Aniline | MIC (µg/mL) | MM/GBSA (kcal/mol) |
|---|---|---|---|
| 20sa23 | 3,5-dichloro | 1 | −32.02 |
| 25sa23 | 3-CF₃, 5-F | 2 | −28.53 |
| 22sa23 | 3-CF₃, 4-F | 4 | −23.35 |
| 5sa23 | 3-fluoro | 16 | −17.89 |
| 17sa23 | 3,4-difluoro | 32 | −10.78 |

MIC values are from the source paper, measured against *Staphylococcus aureus*
ATCC 29213. The structure is *C. difficile* FabK, so the comparison is between
a binding model for one enzyme and whole-cell activity in another organism.
MM/GBSA values are computational scores, not experimentally measured binding
free energies. Five compounds is a small set, and the rank agreement carries
limited statistical weight on its own.

Reference paper: https://doi.org/10.1021/acs.jmedchem.3c00969

---

Siam Chowdhury · www.siamchowdhury.com
Computational and Medicinal Chemistry
Mentor: Mohammad Alam, PhD, Professor of Chemistry
Dr. Alam's Research Team · www.alamresearch.org
Arkansas State University, AR, USA
