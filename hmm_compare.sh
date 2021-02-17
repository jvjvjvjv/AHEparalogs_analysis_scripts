for x in *.fas.fas; do
        header=$(echo $x | sed 's/.fas.fas//')
        HmmCleaner.pl $x -costs -0.75 -0.50 0.15 0.45
        rename .fas_hmm _hmmpcustom2 ./"$header".*
done
