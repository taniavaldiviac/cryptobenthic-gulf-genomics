#!/usr/bin/env bash
# Usage: bash run_batch.sh batch.tsv   (columns: species, run, seed_basename, strategy)
# strategy: WGS (first 5M reads per file), WGSFULL, UCE / exon (full run), ddRAD (mapping mode)
# Seeds must exist as seeds/<seed_basename>.fasta
set -uo pipefail
TSV=${1:?batch.tsv}; LOG=${TSV%.tsv}.log
while IFS=$'\t' read -r sp run seed st; do
  tag=$(echo "$sp" | tr ' ' '_')
  ls asm/${tag}_${run}/*path_sequence.fasta >/dev/null 2>&1 && continue
  echo "== $(date +%T) START $sp $run" >> "$LOG"
  ./run_one.sh "$sp" "$run" "$seed" "$st" >> "$LOG" 2>&1
  [[ "$st" == WGS* || "$st" == ddRAD ]] && rm -f reads/${run}_*
  echo "== $(date +%T) END $sp" >> "$LOG"
done < "$TSV"
echo "== BATCH DONE" >> "$LOG"
