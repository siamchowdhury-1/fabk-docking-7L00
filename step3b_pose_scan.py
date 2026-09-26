#!/usr/bin/env python3
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# STEP 3B: Anchor-contact scan across ALL docked poses.
#
# Pose 1 alone can be luck. This scans every returned pose of
# every compound and asks how often each key FabK anchor is
# engaged. A compound that never reaches an anchor in any pose
# is a far safer call than one judged on a single pose.
#
# Usage:
#   python3 step3b_pose_scan.py                  # global run
#   python3 step3b_pose_scan.py --run focused    # focused run
#   python3 step3b_pose_scan.py --run focused --compounds 5sa23 17sa23
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

# FabK anchors. His143 is catalytic; FMN is the cofactor.
# Side-chain-only atoms are listed where a backbone contact
# would be uninformative (every residue has N, CA, C, O).
ANCHORS = [
    ("FMN",  "FMN 401", None),
    ("HIS143", "HIS 143", {"ND1", "NE2", "CD2", "CE1", "CG"}),
    ("GLU136", "GLU 136", {"OE1", "OE2", "CD"}),
    ("ASN44",  "ASN 44",  {"OD1", "ND2", "CG"}),
    ("MET275", "MET 275", None),
    ("TRP21",  "TRP 21",  None),
]

POLAR_CUT = 3.5     # N/O to N/O, hydrogen-bond-like
CONTACT_CUT = 4.0   # any heavy atom


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run", default="global")
    p.add_argument("--tag", default="7L00_AB")
    p.add_argument("--compounds", nargs="*")
    return p.parse_args()


def read_receptor(path: Path):
    atoms = []
    for l in path.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name = l[12:16].strip()
        el = (l[76:78].strip() if len(l) >= 78 else "") or re.sub(r"[^A-Za-z]", "", name)[:1]
        if el.upper() == "H":
            continue
        atoms.append({
            "name": name,
            "res": l[17:20].strip(),
            "chain": l[21],
            "num": l[22:26].strip(),
            "el": el.upper(),
            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
        })
    return atoms


def read_poses(sdf: Path, pdbqt: Path):
    blocks, cur = [], []
    for line in sdf.read_text().splitlines():
        if line.startswith("$$$$"):
            blocks.append(cur)
            cur = []
        else:
            cur.append(line)

    poses = []
    for rec in blocks:
        if len(rec) < 4:
            continue
        n = int(rec[3][:3])
        atoms = [
            {"el": l[31:34].strip(),
             "xyz": (float(l[0:10]), float(l[10:20]), float(l[20:30]))}
            for l in rec[4:4 + n]
        ]
        poses.append({"atoms": [a for a in atoms if a["el"] != "H"], "score": None})

    scores = [
        float(l.split()[3]) for l in pdbqt.read_text().splitlines()
        if l.startswith("REMARK VINA RESULT:")
    ]
    for i, p in enumerate(poses):
        if i < len(scores):
            p["score"] = scores[i]
    return poses


def anchor_state(pose, rec_by_anchor):
    """Return per-anchor (contact, polar_pairs) for one pose."""
    out = {}
    for label, atoms in rec_by_anchor.items():
        contact, polar = False, 0
        for r in atoms:
            for a in pose["atoms"]:
                d = math.dist(r["xyz"], a["xyz"])
                if d <= CONTACT_CUT:
                    contact = True
                if d <= POLAR_CUT and r["el"] in "NO" and a["el"] in "NO":
                    polar += 1
        out[label] = (contact, polar)
    return out


def main():
    args = parse_args()
    run = args.run
    prefix = "docking_" if run == "global" else f"docking_{run}_"
    txt_dir = Path("txt_outputs") if run == "global" else Path("txt_outputs") / run
    label = "" if run == "global" else f" [{run}]"

    rec_path = Path(f"{args.tag}_protein_FMN.pdb")
    if not rec_path.is_file():
        raise SystemExit(f"ERROR: {rec_path} not found. Run Step 2 first.")

    cpds = args.compounds
    if not cpds:
        cpds = sorted(
            d.name[len(prefix):] for d in Path().glob(f"{prefix}*") if d.is_dir()
        )
    if not cpds:
        raise SystemExit(f"ERROR: no {prefix}* folders found.")

    rec = read_receptor(rec_path)
    rec_by_anchor = {}
    for label_a, spec, side in ANCHORS:
        res, num = spec.split()
        sel = [
            a for a in rec
            if a["res"] == res and a["num"] == num
            and (side is None or a["name"] in side)
        ]
        if sel:
            rec_by_anchor[label_a] = sel

    missing = [a[0] for a in ANCHORS if a[0] not in rec_by_anchor]

    head = ["=" * 68, *IDENTITY,
            f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}", "=" * 68, ""]

    lines = head + [
        f"34 ANCHOR-CONTACT SCAN ACROSS ALL POSES{label}",
        "-" * 68,
        f"Receptor : {rec_path}",
        f"Contact  : any heavy atom within {CONTACT_CUT:.1f} A",
        f"Polar    : ligand N/O within {POLAR_CUT:.1f} A of anchor N/O",
        "His143, Glu136 and Asn44 are counted on SIDE-CHAIN atoms only,",
        "so a backbone brush does not count as engaging the residue.",
        "",
    ]
    if missing:
        lines += [f"NOTE: anchors not found in receptor: {', '.join(missing)}", ""]

    anchors = list(rec_by_anchor)
    summary = {}

    for c in cpds:
        sdf = Path(f"{prefix}{c}/{c}_poses.sdf")
        qt = Path(f"{prefix}{c}/{c}_out.pdbqt")
        if not (sdf.is_file() and qt.is_file()):
            lines += [f"{c}: pose files missing, skipped", ""]
            continue

        poses = read_poses(sdf, qt)
        states = [anchor_state(p, rec_by_anchor) for p in poses]

        lines += [
            f"COMPOUND {c}   ({len(poses)} poses)",
            "-" * 68,
            f"{'POSE':<6}{'SCORE':>8}  " + "".join(f"{a:>10}" for a in anchors),
        ]
        for i, (p, st) in enumerate(zip(poses, states), start=1):
            cells = []
            for a in anchors:
                contact, polar = st[a]
                cells.append(f"{('%d H' % polar) if polar else ('yes' if contact else '-'):>10}")
            lines.append(f"{i:<6}{p['score']:>8.2f}  " + "".join(cells))

        counts = {
            a: (
                sum(1 for st in states if st[a][0]),
                sum(1 for st in states if st[a][1] > 0),
            )
            for a in anchors
        }
        summary[c] = (len(poses), counts, poses[0]["score"])

        lines += [
            "",
            f"{'poses with contact':<22}" + "".join(f"{counts[a][0]:>10}" for a in anchors),
            f"{'poses with H-bond':<22}" + "".join(f"{counts[a][1]:>10}" for a in anchors),
            "",
        ]

    if summary:
        lines += [
            "=" * 68,
            "SUMMARY: how many poses engage each anchor (contact / H-bond)",
            "=" * 68,
            f"{'COMPOUND':<11}{'POSES':>6}{'BEST':>8}" + "".join(f"{a:>12}" for a in anchors),
        ]
        for c, (n, counts, best) in summary.items():
            lines.append(
                f"{c:<11}{n:>6}{best:>8.2f}"
                + "".join(f"{counts[a][0]}/{counts[a][1]}".rjust(12) for a in anchors)
            )

        never = {
            c: [a for a in anchors if counts[a][0] == 0]
            for c, (_, counts, _) in summary.items()
        }
        lines += ["", "ANCHORS NEVER REACHED IN ANY POSE", "-" * 68]
        for c, miss in never.items():
            lines.append(f"  {c:<11}{', '.join(miss) if miss else 'none - reaches every anchor'}")

        lines += [
            "",
            "READING THIS TABLE",
            "-" * 68,
            "A compound that reaches an anchor in only its best pose may have",
            "been lucky. One that reaches it across many poses is engaging that",
            "residue robustly. One that never reaches an anchor in 20 tries is",
            "the strongest negative signal rigid docking can give.",
            "",
            "Vina scores are computational docking scores,",
            "not experimentally measured binding free energies.",
        ]

    txt_dir.mkdir(parents=True, exist_ok=True)
    out = txt_dir / "34_anchor_contact_pose_scan.txt"
    out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines[len(head):]))
    print(f"\nWritten: {out}")


if __name__ == "__main__":
    main()
