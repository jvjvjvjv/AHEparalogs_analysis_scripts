#!/bin/bash
#SBATCH --job-name=blastx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=30G
#SBATCH --mail-user=jason.vailionis@uconn.edu
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err

echo $(hostname)
date

module load blast

laopath="/home/FCAM/egordon/genomes/Laodelphax/Laoproteome.fasta"
nilapath="/home/FCAM/egordon/genomes/AHEloci/Nilproteome.fasta"

mkdir to_laoprot
mkdir to_nilaprot

for x in *_subset.fasta; do
        blastx -query $x -db $laopath -num_threads 20 -outfmt 6 -out ./to_laoprot/"$x".res
        sort -u -k1,1 ./to_laoprot/"$x".res > ./to_laoprot/"$x"_tophits.res
done

for x in *_subset.fasta; do
        blastx -query $x -db $nilapath -num_threads 20 -outfmt 6 -out ./to_nilaprot/"$x".res
        sort -u -k1,1 ./to_nilaprot/"$x".res > ./to_nilaprot/"$x"_tophits.res
done
