#!/usr/bin/env bash
# Recupera mitogenomas de lecturas publicas de SRA para criptobentonicos del Golfo de California.
# Uso: bash sra_mitogenomas.sh corridas_sra.tsv [hilos]
# Modo "asm": ensamblaje con GetOrganelle (WGS, UCE, captura de exones).
# Modo "map": mapeo a semilla del mismo genero y consenso (ddRAD; cobertura parcial esperada).
set -euo pipefail
TSV=${1:-corridas_sra.tsv}; T=${2:-6}; WGS_SPOTS=5000000
mkdir -p seeds reads asm map logs
tail -n +2 "$TSV" | while IFS=$'\t' read -r sp fam run strat modo seed gb ref nota; do
  tag=$(echo "$sp" | tr ' ' '_'); log=logs/${tag}_${run}.log
  [ -s seeds/$seed.fasta ] || efetch -db nuccore -id "$seed" -format fasta > seeds/$seed.fasta
  if [ ! -s reads/${run}_1.fastq.gz ]; then
    if [ "$strat" = "WGS" ]; then X="-X $WGS_SPOTS"; else X=""; fi
    fastq-dump --split-files --gzip $X -O reads "$run" >>"$log" 2>&1 || { echo "falla descarga $run"; continue; }
  fi
  R1=reads/${run}_1.fastq.gz; R2=reads/${run}_2.fastq.gz
  if [ "$modo" = "asm" ]; then
    get_organelle_from_reads.py -1 "$R1" -2 "$R2" -s seeds/$seed.fasta -F animal_mt \
      -R 10 -k 21,45,65,85 -t "$T" -o asm/${tag}_${run} --overwrite >>"$log" 2>&1 \
      || echo "GetOrganelle sin grafo completo: $run (revisar; usar modo map como respaldo)"
  else
    [ -s seeds/$seed.fasta.bwt ] || bwa index seeds/$seed.fasta >/dev/null 2>&1
    bwa mem -t "$T" seeds/$seed.fasta "$R1" "$R2" 2>>"$log" \
      | samtools view -b -F 4 -q 20 - | samtools sort -o map/${tag}_${run}.bam -
    samtools index map/${tag}_${run}.bam
  fi
done
# Consenso por especie en modo map (fusiona individuos de la misma especie)
for sp in $(ls map/*.bam 2>/dev/null | sed 's#map/##; s#_SRR.*##' | sort -u); do
  samtools merge -f map/${sp}_merged.bam map/${sp}_SRR*.bam && samtools index map/${sp}_merged.bam
  samtools consensus -a --show-ins no --min-depth 3 -f fasta map/${sp}_merged.bam -o map/${sp}_consenso.fasta
  samtools depth -a map/${sp}_merged.bam | awk -v s=$sp '{n++; if($3>=3)c++} END{printf "%s\tposiciones=%d\tcubiertas_3x=%.1f%%\n", s, n, 100*c/n}' >> map/resumen_cobertura.tsv
done
echo "Listo. Ensamblajes en asm/, consensos en map/, cobertura en map/resumen_cobertura.tsv"
