#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=ssp_seq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/seq_%j.out
#SBATCH --error=logs/seq_%j.err

set -euo pipefail

BIN=./shsup_seq
CSV=results/seq_times.csv
REPS=${REPS:-20}

mkdir -p logs results

# cabeçalho do CSV, se ainda não existir
if [ ! -s "$CSV" ]; then
  echo "input,tempo,sd" > "$CSV"
fi

# percorre todos os arquivos em ./inputs
for INPUT in $(ls inputs/*.txt | sort -V); do
  INBASE=$(basename "$INPUT")

  TIMES=()
  for ((r=1; r<=REPS; r++)); do
    # programa imprime superstring + tempo; pegamos só a última linha
    T=$("$BIN" < "$INPUT" | tail -n 1)
    TIMES+=("$T")
  done

  printf '%s\n' "${TIMES[@]}" | \
  awk -v input="$INBASE" '
    { x[NR]=$1; s+=$1 }
    END {
      n=NR
      m=s/n
      for (i=1;i<=n;i++) { d=x[i]-m; v+=d*d }
      sd=(n>1)?sqrt(v/(n-1)):0
      printf "%s,%f,%f\n", input, m, sd
    }
  ' >> "$CSV"
done
