awk -v FS="\t" -v OFS="\t" '{print($2, ($9-1), $10, $1, $11)}' I6843contigs_subset.fasta_tophits.res | sort-bed - > I6843contigs.bed
