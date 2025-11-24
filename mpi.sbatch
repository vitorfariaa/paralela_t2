#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=ssp_mpi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=6
#SBATCH --time=02:00:00
#SBATCH --output=logs/mpi_%j.out
#SBATCH --error=logs/mpi_%j.err

set -euo pipefail

BIN=./shsup_mpi
CSV=results/mpi_times.csv
REPS=${REPS:-20}

# lista de np que você quer testar
NP_LIST=(1 2 3 4 6 8 10 12)

mkdir -p logs results

if [ ! -s "$CSV" ]; then
  echo "input,np,tempo,sd" > "$CSV"
fi

for INPUT in $(ls inputs/*.txt | sort -V); do
  INBASE=$(basename "$INPUT")

  for NP in "${NP_LIST[@]}"; do
    # processos por nó (ceil(NP/nnodes))
    PPR=$(( (NP + SLURM_NNODES - 1) / SLURM_NNODES ))

    TIMES=()
    for ((r=1; r<=REPS; r++)); do
      # srun chama o binário MPI; pegamos só o tempo (última linha)
      T=$(srun --ntasks="$NP" --ntasks-per-node="$PPR" "$BIN" < "$INPUT" | tail -n 1)
      TIMES+=("$T")
    done

    printf '%s\n' "${TIMES[@]}" | \
    awk -v input="$INBASE" -v np="$NP" '
      { x[NR]=$1; s+=$1 }
      END {
        n=NR
        m=s/n
        for (i=1;i<=n;i++) { d=x[i]-m; v+=d*d }
        sd=(n>1)?sqrt(v/(n-1)):0
        printf "%s,%d,%f,%f\n", input, np, m, sd
      }
    ' >> "$CSV"
  done
done
