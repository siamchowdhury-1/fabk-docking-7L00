#!/usr/bin/env bash
# ============================================================
# Siam Chowdhury (www.siamchowdhury.com)
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Arkansas State University, AR, USA
# ------------------------------------------------------------
# SHOW ALL RESULTS
#
# Gathers every .txt report and .png graph already produced and
# builds ONE self-contained HTML file. Images are embedded, so
# the file opens anywhere with nothing else attached.
#
# Nothing is moved, renamed or deleted. Originals stay put.
#
#   ./show_all_results.sh              build the report and open it
#   ./show_all_results.sh --print      also dump every txt to the terminal
#   ./show_all_results.sh --no-open    build only
# ============================================================

set -Eeuo pipefail

OUT="7L00_FabK_results_report.html"
PRINT=0
OPEN=1

for arg in "$@"; do
    case "$arg" in
        --print)   PRINT=1 ;;
        --no-open) OPEN=0 ;;
        *) echo "unknown option: $arg"; exit 1 ;;
    esac
done

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found."; exit 1; }

echo "Collecting results in $(pwd) ..."

python3 - "$OUT" <<'PY'
import base64
import html
import re
import sys
from datetime import datetime
from pathlib import Path

out_path = Path(sys.argv[1])
root = Path(".")

IDENTITY = [
    "siam chowdhury",
    "Computational and Medicinal Chemistry",
    "[Dr. Alam's Research Team] www.alamresearch.org",
    "Arkansas State University",
]

# Section order. Anything not matched lands in "Other files".
SECTIONS = [
    ("Structure preparation", ["01_", "02_", "03_", "04_", "05_", "06_", "07_"]),
    ("Binding site and receptor", ["08_", "09_", "10_", "11_", "12_", "14_"]),
    ("Redocking validation (XCJ)", ["13_", "15_", "16_", "17_", "17b_", "18_"]),
    ("Compound preparation", ["19_", "20_", "21_"]),
    ("Docking runs", ["22_", "23_", "24_"]),
    ("Contacts and ranking", ["25_", "26_", "27_", "34_"]),
    ("MM/GBSA", ["35_", "36_", "37_", "38_"]),
    ("Terminal logs", ["00_"]),
]


def run_label(path: Path) -> str:
    parts = path.parts
    if "focused" in parts:
        return "focused (16 A box)"
    if parts[0] == "txt_outputs":
        return "global (22 A box)"
    return ""


txts = sorted(
    p for p in root.rglob("*.txt")
    if p.is_file()
    and "txt_outputs" in p.parts
    and p.stat().st_size > 0
)
pngs = sorted(
    p for p in root.rglob("*.png")
    if p.is_file() and p.parts[0].startswith("graphs")
)

if not txts and not pngs:
    raise SystemExit("ERROR: no txt_outputs/*.txt or graphs*/*.png files found. "
                     "Run this from your working folder.")


def section_of(path: Path) -> str:
    name = path.name
    for title, prefixes in SECTIONS:
        if any(name.startswith(p) for p in prefixes):
            return title
    return "Other files"


def strip_identity(text: str) -> str:
    """Remove the repeated identity header; it is shown once at the top."""
    lines = text.splitlines()
    rule = lambda l: set(l.strip()) == {"="} and len(l.strip()) > 20

    # The header is: rule, identity lines, rule. Drop the whole block.
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

NOW = f"{datetime.now():%d %B %Y, %H:%M}"

parts = ["""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>7L00 FabK docking results</title>
<style>
  :root { --ink:#1a1a1a; --muted:#666; --line:#d8d8d8; --bg:#fff; --card:#fafafa; }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--bg); color:var(--ink);
         font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
  .wrap { max-width:1100px; margin:0 auto; padding:0 24px 80px; }
  header { border-bottom:3px solid var(--ink); padding:32px 0 20px; margin-bottom:8px; }
  h1 { margin:0 0 6px; font-size:26px; letter-spacing:-.01em; }
  .who { color:var(--muted); font-size:14px; line-height:1.5; }
  .who strong { color:var(--ink); }
  nav { position:sticky; top:0; background:var(--bg); padding:14px 0;
        border-bottom:1px solid var(--line); margin-bottom:28px; z-index:10; }
  nav a { display:inline-block; margin:3px 14px 3px 0; color:#0b5; text-decoration:none;
          font-size:13px; font-weight:600; }
  nav a:hover { text-decoration:underline; }
  h2 { margin:44px 0 4px; font-size:20px; padding-bottom:8px;
       border-bottom:2px solid var(--line); }
  details { border:1px solid var(--line); border-radius:6px; margin:10px 0;
            background:var(--card); }
  details[open] { background:#fff; }
  summary { cursor:pointer; padding:11px 14px; font-weight:600; font-size:14px;
            list-style:none; display:flex; justify-content:space-between; gap:12px; }
  summary::-webkit-details-marker { display:none; }
  summary::before { content:"\\25B8"; margin-right:9px; color:var(--muted);
                    transition:transform .15s; }
  details[open] summary::before { content:"\\25BE"; }
  .path { color:var(--muted); font-weight:400; font-size:12px;
          font-family:ui-monospace,Menlo,Consolas,monospace; }
  pre { margin:0; padding:14px 16px; overflow-x:auto; font-size:12.5px; line-height:1.5;
        font-family:ui-monospace,Menlo,Consolas,monospace; border-top:1px solid var(--line);
        background:#fff; white-space:pre; }
  figure { margin:0 0 26px; }
  figure img { width:100%; border:1px solid var(--line); border-radius:6px; }
  figcaption { color:var(--muted); font-size:13px; margin-top:7px; }
  .grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(420px,1fr)); gap:24px; }
  .tag { display:inline-block; background:#eef; color:#335; border-radius:4px;
         padding:1px 7px; font-size:11px; font-weight:600; margin-left:8px; }
  .bar { display:flex; gap:10px; margin:18px 0 6px; flex-wrap:wrap; }
  button { font:inherit; font-size:13px; padding:7px 14px; border:1px solid var(--line);
           background:#fff; border-radius:6px; cursor:pointer; }
  button:hover { background:var(--card); }
  footer { margin-top:56px; padding-top:18px; border-top:1px solid var(--line);
           color:var(--muted); font-size:12.5px; }
  @media print {
    nav, .bar { display:none; }
    details { break-inside:avoid; }
    details:not([open]) > *:not(summary) { display:block; }
  }
</style></head><body><div class="wrap">
<header>
<h1>Molecular docking: 7L00 <em>C. difficile</em> FabK</h1>
<div class="who">"""
+ "<br>".join(
    f"<strong>{html.escape(IDENTITY[0])}</strong>" if i == 0 else html.escape(l)
    for i, l in enumerate(IDENTITY)
)
+ f"""<br>Report generated {NOW}</div>
</header>
<nav>"""]

nav = []
for title, _ in SECTIONS + [("Other files", [])]:
    if title in grouped:
        nav.append(f'<a href="#{re.sub(r"[^a-z]+", "-", title.lower())}">{html.escape(title)}</a>')
if pngs:
    nav.insert(0, '<a href="#graphs">Graphs</a>')
parts.append("".join(nav))
parts.append("""</nav>
<div class="bar">
  <button onclick="document.querySelectorAll('details').forEach(d=>d.open=true)">Expand all</button>
  <button onclick="document.querySelectorAll('details').forEach(d=>d.open=false)">Collapse all</button>
  <button onclick="window.print()">Print / save as PDF</button>
</div>""")

# ---- graphs first: they carry a meeting ----
if pngs:
    parts.append('<h2 id="graphs">Graphs</h2>')
    for folder in sorted(png_groups):
        tag = "focused (16 A box)" if "focused" in folder else "global (22 A box)"
        parts.append(f'<h3 style="font-size:16px;color:#666;margin:26px 0 12px">'
                     f'{html.escape(folder)}/ <span class="tag">{tag}</span></h3>')
        parts.append('<div class="grid">')
        for p in png_groups[folder]:
            b64 = base64.b64encode(p.read_bytes()).decode()
            parts.append(
                f'<figure><img src="data:image/png;base64,{b64}" alt="{html.escape(p.name)}">'
                f'<figcaption>{html.escape(p.name)}</figcaption></figure>'
            )
        parts.append("</div>")

# ---- text reports ----
for title, _ in SECTIONS + [("Other files", [])]:
    if title not in grouped:
        continue
    anchor = re.sub(r"[^a-z]+", "-", title.lower())
    parts.append(f'<h2 id="{anchor}">{html.escape(title)}</h2>')
    for p in sorted(grouped[title]):
        body = strip_identity(p.read_text(errors="replace"))
        label = run_label(p)
        tag = f'<span class="tag">{html.escape(label)}</span>' if label else ""
        parts.append(
            "<details><summary><span>" + html.escape(p.name) + tag +
            f'</span><span class="path">{html.escape(str(p))}</span></summary>'
            f"<pre>{html.escape(body)}</pre></details>"
        )

parts.append(
    "<footer>" + " &middot; ".join(html.escape(l) for l in IDENTITY) +
    f"<br>{len(txts)} text reports and {len(pngs)} graphs, "
    "collected from their original locations; no file was moved or altered."
    "<br>Vina and MM/GBSA values are computational scores, "
    "not experimentally measured binding free energies.</footer>"
    "</div></body></html>"
)

out_path.write_text("".join(parts), encoding="utf-8")
size = out_path.stat().st_size / 1024
print(f"{len(txts)} text reports, {len(pngs)} graphs -> {out_path} ({size:.0f} KB)")
PY

if (( PRINT )); then
    echo
    echo "############################################################"
    echo "# ALL TEXT REPORTS"
    echo "############################################################"
    find txt_outputs -name '*.txt' -size +0 | sort | while read -r f; do
        echo
        echo "============================================================"
        echo "FILE: $f"
        echo "============================================================"
        cat "$f"
    done
fi

echo
echo "Report ready: $(pwd)/$OUT"

if (( OPEN )); then
    if command -v wslview >/dev/null 2>&1; then
        wslview "$OUT" 2>/dev/null &
    elif command -v explorer.exe >/dev/null 2>&1; then
        explorer.exe "$(wslpath -w "$OUT")" 2>/dev/null &
    elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open "$OUT" 2>/dev/null &
    else
        echo "Open it manually in your browser."
    fi
    sleep 1
fi
