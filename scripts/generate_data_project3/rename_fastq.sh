sed 's/,/ /g' lookup_bc_organism.csv | while read BC NAME
do
    mv "$BC".fastq.gz "$NAME.fastq.gz"
done