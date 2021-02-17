xargs -a test_batch -i cp {} .  #cp alex's files here
bash hmm_compare.sh
mkdir cleaned_fastas
mv *.fasta cleaned_fastas/
cd cleaned_fastas/
Rscript /home/CAM/jvailionis/scripts/view_DNAalignments.R -d . 12  #mine is edited to not have a custom library location
mv loci1.png hmmtest.png
