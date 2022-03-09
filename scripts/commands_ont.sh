
cd ~/workdir/groupwork_ont/reads
BASENAMES=`ls *.fastq.gz | cut -f 1 -d "."`

cd ..

minimap2 \
-x splice \
-d reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa

for name in $BASENAMES
do
  minimap2 \
  -a \
  -x splice \
  -G 500k \
  -t 4 \
  reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
  reads/$name.fastq.gz \
  | samtools sort \
  | samtools view -bh > alignments/$name.bam
done

samtools merge alignments/merged.bam alignments/*.bam
samtools index alignments/merged.bam

git clone https://github.com/BrooksLabUCSC/flair.git

python flair/bin/bam2Bed12.py \
-i alignments/merged.bam \
> alignments/merged.bed

mkdir flair_output

python flair/flair.py correct \
-q alignments/merged.bed \
-g reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-f reference/Homo_sapiens.GRCh38.102.gtf \
-o flair_output/CACNA1C

READS=`ls reads/*.fastq.gz | tr "\n" ","`
READS="${READS%?}" #remove last comma

python flair/flair.py collapse \
-g reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-f reference/Homo_sapiens.GRCh38.102.gtf \
-r $READS \
-q flair_output/CACNA1C_all_corrected.bed \
-o flair_output/CACNA1C.collapse

python flair/flair.py quantify \
-r reads_manifest.tsv \
-i flair_output/CACNA1C.collapse.isoforms.fa

python flair/bin/plot_isoform_usage.py \
flair_output/CACNA1C.collapse.isoforms.bed \
counts_matrix.tsv \
ENSG00000151067
