#!/bin/bash
#SBATCH --job-name=blast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 3
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=5G
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err

module load blast

for x in /home/FCAM/egordon/genomes/AHEloci/Heteropterapslit/processed/test/*.fasta.fas; do
        for y in /home/CAM/mstukel/AHE/Kikihia/blastdb/*.fasta; do
                tblastn -query $x -db $y -max_target_seqs 5 -outfmt 6 -num_threads 3 | cut -f 11 | sed "s|^|$x\t$y\t|" | awk -F '\t' -v OFS='\t' '{sub(/.*\//, "", $1); sub(/.*\//, "", $2)} 1' >> 2_15_21_blastscores_test.tsv
        done
done
