# script to split a fasta file into a ton of files and return the file paths so that you can pipe it to parallel
inputfile=`readlink -f $1`
outdir="${1%.*}"
mkdir $outdir
cd $outdir

while read line ; do
  if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < $inputfile

ls "$PWD/"*
