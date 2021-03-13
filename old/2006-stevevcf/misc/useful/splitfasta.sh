# dc into folder of orig filename
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        #outfile=${line#>}.fa
	outfile=${line//\//\_}
	outfile=${outfile#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < *.fasta



${parameter/"\?(.*)"/""}:

