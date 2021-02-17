awk -v FS="\t" -v OFS="\t" '{print($2, ($9-1), $10, $1, $11)}' I6843contigs_subset.fasta_tophits.res | sort-bed - > I6843contigs.bed

awk '$5 < 1e-10' I6843contigs.bed | bedtools merge #get full intervals from overlap regions

awk '$5 < 1e-10' I6843contigs.bed | bedtools cluster | cut -f 6 | uniq -c # get counts of overlap clusters

for file in *tophits.res; do
  outname="$(echo $file | sed 's/_subset.fasta_tophits.res//')"
  awk -v FS="\t" -v OFS="\t" '{print($2, ($9-1), $10, $1, $11)}' $file | awk '$5 < 1e-20' | sort-bed - > "$outname".bed
done

#lol this wont work
# bedtools intersect -a "$(ls *.bed | head -n 1)" -b "$(ls *.bed | tail -n +2 | tr '\n' ',' | sed 's/,$//' | sed 's/,/, /g')"


cp /home/FCAM/egordon/genomes/Laodelphax/Laoproteome.fasta .

#get all the regions with sequences from at least 4 out of 5 assemblies, then gets the sequence, min length 30
multiIntersectBed -i *.bed | awk '$4 > 3' | bedtools getfasta -fi Laoproteome.fasta -bed stdin | seqkit seq -m 30 > laoprothits.fasta
#the intersect gets the intersecting region so you get small bits, i want to merge instead

cat *contigs.bed | awk '$5 < 1e-10' | sort-bed - | bedtools cluster | cut -f 6 | uniq -c

cat *contigs.bed | awk '$5 < 1e-10' | sort-bed - | bedtools merge | bedtools coverage -a stdin -b <(cat *contigs.bed | awk '$5 < 1e-10') | cut -f 1,2,3,4

tblastn -query laoprothits.fasta -db /home/FCAM/egordon/genomes/Laodelphax/Laogenome.fasta -outfmt 6 -max_target_seqs 1 | awk -v FS="\t" -v OFS="\t" '{print($2, ($9-1), $10, $1, $11)}' > laogenomehits.bed

#why are so many of the proteome to genome hits weak? were they just short queries?

#original number of sequences in 5 subsets = 10677
seqkit stats ../*subset.fasta | awk '{print $4}' | tail -n +2 | paste -sd+ | tr -d ',' | bc

cat *contigs.bed | sort-bed - | bedtools merge | wc

wc *contigs.bed
