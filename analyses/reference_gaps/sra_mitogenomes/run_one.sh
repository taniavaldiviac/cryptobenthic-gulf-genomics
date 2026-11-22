#!/usr/bin/env bash
# run1.sh ESPECIE RUN SEMILLA ESTRATEGIA  (descarga de ENA + GetOrganelle o mapeo)
set -uo pipefail
[ -n "${SRAMITO_ENV:-}" ] && export PATH=$SRAMITO_ENV/bin:$PATH
sp=$1; run=$2; seed=$3; strat=$4; tag=$(echo "$sp" | tr ' ' '_')
mkdir -p reads asm map logs
urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2 | tr ';' ' ')
for u in $urls; do b=$(basename $u); case $b in *_1.fastq.gz|*_2.fastq.gz) ;; *) continue;; esac; f=reads/$b
  if [ ! -s $f ]; then
    if [ "$strat" = WGS ]; then curl -s "https://$u" | gzip -dc 2>/dev/null | head -n 20000000 | gzip -1 > $f
    else curl -s -o $f "https://$u"; fi
  fi; done
R1=reads/${run}_1.fastq.gz; R2=reads/${run}_2.fastq.gz
n1=$(gzip -dc $R1 2>/dev/null | wc -l); n2=$(gzip -dc $R2 2>/dev/null | wc -l); n=$(( (n1<n2?n1:n2)/4*4 ))
if [ $n1 -ne $n2 ]; then echo "pares desiguales $n1/$n2 -> recorto a $n"; for f in $R1 $R2; do gzip -dc $f 2>/dev/null | head -n $n | gzip -1 > $f.tmp && mv $f.tmp $f; done; fi
echo "$run lecturas: $(( $(gzip -dc $R1 | wc -l) / 4 ))"
if [ "$strat" = ddRAD ]; then
  [ -s seeds/$seed.fasta.bwt ] || bwa index seeds/$seed.fasta 2>/dev/null
  bwa mem -t 2 seeds/$seed.fasta $R1 $R2 2>/dev/null | samtools view -b -F 4 -q 20 - | samtools sort -o map/${tag}_${run}.bam - && samtools index map/${tag}_${run}.bam
  echo "mapeadas: $(samtools view -c map/${tag}_${run}.bam)"
else
  get_organelle_from_reads.py -1 $R1 -2 $R2 -s seeds/$seed.fasta -F animal_mt -R 15 -t 2 --memory-save --max-n-words 2E7 -J 3 -M 5 -k 21,55,85 --spades-options="--only-assembler" --reduce-reads-for-coverage inf --max-reads inf \
    -o asm/${tag}_${run} --continue > logs/${tag}_${run}.log 2>&1
  tail -4 logs/${tag}_${run}.log; ls asm/${tag}_${run}/*.fasta 2>/dev/null
fi
