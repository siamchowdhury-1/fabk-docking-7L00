#!/usr/bin/env bash
# Siam Chowdhury | [Dr. Alam's Research Team] | Arkansas State University
# Redraw the COMBINED graphs for the five-compound series.
# The compound name survives in these PNGs as pixels, so no text
# search can see it. GRAPHS_ONLY=1 means no docking and no ligand
# preparation, so the scores cannot move - and this script proves
# it by fingerprinting every score file before and after.
set -Eeuo pipefail
RUN="${RUN:-focused}"; BOX="${BOX:-16}"; APPLY=0
for a in "$@"; do case "$a" in --apply) APPLY=1;; *) echo "unknown: $a"; exit 1;; esac; done
COMPOUNDS=(20sa23 25sa23 22sa23 5sa23 17sa23)
GD=$([ "$RUN" = global ] && echo graphs || echo "graphs_${RUN}")
PP=$([ "$RUN" = global ] && echo docking_ || echo "docking_${RUN}_")
TD=$([ "$RUN" = global ] && echo txt_outputs || echo "txt_outputs/${RUN}")
COMBINED=("$GD/multi_docking_7L00_FabK.png" "$GD/best_scores_7L00_FabK.png" "$GD/MIC_vs_score_7L00_FabK.png")
echo "============================================================"
echo "siam chowdhury | [Dr. Alam's Research Team]"
echo "Arkansas State University | $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo; echo "REDRAW COMBINED GRAPHS: ${COMPOUNDS[*]}"
[ "$APPLY" = 1 ] && echo "MODE: APPLY" || echo "MODE: DRY RUN (add --apply)"
echo
WATCH=()
for c in "${COMPOUNDS[@]}"; do
  for f in "$PP$c/${c}_summary.tsv" "$PP$c/${c}_out.pdbqt" "$PP$c/${c}.pdbqt" \
           "$PP$c/${c}_poses.sdf" "$TD/23_${c}_scores.txt" "$TD/22_${c}_vina_run.txt"; do
    [ -f "$f" ] && WATCH+=("$f")
  done
done
[ -f mic_values.txt ] && WATCH+=(mic_values.txt)
[ ${#WATCH[@]} -eq 0 ] && { echo "ERROR: no score files found. Wrong directory?"; exit 1; }
BEFORE=$(mktemp); for f in "${WATCH[@]}"; do md5sum "$f"; done | sort -k2 > "$BEFORE"
echo "PROTECTED FILES: ${#WATCH[@]} fingerprinted"; echo
echo "SCORES BEFORE"; echo "------------------------------------------------------------"
for c in "${COMPOUNDS[@]}"; do
  s="$PP$c/${c}_summary.tsv"
  [ -f "$s" ] && printf '  %-10s %s\n' "$c" "$(awk 'NR>1{print $2; exit}' "$s")" \
              || printf '  %-10s summary missing\n' "$c"
done
echo
echo "COMBINED GRAPHS"; echo "------------------------------------------------------------"
for g in "${COMBINED[@]}"; do
  [ -f "$g" ] && printf '  %-50s %s  %s\n' "$g" "$(du -h "$g"|cut -f1)" "$(date -r "$g" '+%Y-%m-%d %H:%M')" \
              || printf '  %-50s not present\n' "$g"
done
echo
echo "PER-COMPOUND GRAPHS"; echo "------------------------------------------------------------"
for g in "$GD"/Graph_*.png; do
  [ -f "$g" ] || continue
  c=$(basename "$g" .png | sed 's/^Graph_7L00_FabK_//'); k=no
  for m in "${COMPOUNDS[@]}"; do [ "$c" = "$m" ] && k=yes; done
  [ "$k" = yes ] && printf '  %-50s keep\n' "$g" || printf '  %-50s NOT IN THE SERIES\n' "$g"
done
echo
if [ "$APPLY" != 1 ]; then rm -f "$BEFORE"; echo "Dry run complete. Nothing changed."; exit 0; fi
echo "DELETING"; echo "------------------------------------------------------------"
for g in "${COMBINED[@]}"; do [ -f "$g" ] && { rm -f "$g"; echo "  removed  $g"; }; done
echo
echo "REDRAWING (GRAPHS_ONLY=1: no docking)"; echo "------------------------------------------------------------"
LOG=$(mktemp); set +e
GRAPHS_ONLY=1 RUN="$RUN" BOX="$BOX" ./step3_7L00_dock_compounds.sh "${COMPOUNDS[@]}" > "$LOG" 2>&1
RC=$?; set -e
[ $RC -ne 0 ] && { echo "  step3 exited $RC, last lines:"; tail -15 "$LOG" | sed 's/^/    /'; }
echo
echo "VERIFICATION: did any score change?"; echo "------------------------------------------------------------"
AFTER=$(mktemp); for f in "${WATCH[@]}"; do [ -f "$f" ] && md5sum "$f" || echo "MISSING $f"; done | sort -k2 > "$AFTER"
if diff -q "$BEFORE" "$AFTER" >/dev/null 2>&1; then
  echo "  PASS: every score file byte-for-byte identical."
else
  echo "  FAIL: these changed - NOT score-safe:"
  diff "$BEFORE" "$AFTER" | grep '^[<>]' | awk '{print "    " $3}' | sort -u
  echo "  Restore with: git checkout -- $PP* $TD mic_values.txt"
fi
echo
echo "SCORES AFTER"; echo "------------------------------------------------------------"
for c in "${COMPOUNDS[@]}"; do
  s="$PP$c/${c}_summary.tsv"
  [ -f "$s" ] && printf '  %-10s %s\n' "$c" "$(awk 'NR>1{print $2; exit}' "$s")"
done
echo
echo "NEW GRAPHS"; echo "------------------------------------------------------------"
for g in "${COMBINED[@]}"; do
  [ -f "$g" ] && printf '  %-50s %s  %s\n' "$g" "$(du -h "$g"|cut -f1)" "$(date -r "$g" '+%Y-%m-%d %H:%M')" \
              || printf '  %-50s NOT REDRAWN - check log\n' "$g"
done
rm -f "$BEFORE" "$AFTER" "$LOG"
echo
echo "Open each new graph and count the bars: five, not six."
