# lima \
# --ccs \
# --same \
# --split-bam-named \
# --peek-guess \
# m64012_191221_044659.ccsset.bam \
# Barcoded_Adapter_8B.fasta \
# m64012_191221_044659.demux.bam

# NanoPlot \
# --ubam reads/m64012_191221_044659.demux.bc1020--bc1020.bam \
# --outdir nanoplot_reports

mkdir alignments

for sample in `seq 1015 1022`
do
  minimap2 \
  -a \
  -x asm20 \
  reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  reads/$sample.fastq.gz \
  | samtools sort \
  | samtools view -bh \
  > alignments/"$sample".bam

  samtools index alignments/"$sample".bam
done

git clone https://github.com/PacificBiosciences/apps-scripts.git

for gene in gene1 gene2
do
  apps-scripts/RepeatAnalysisTools/makeReports.sh \
  targets/target_"$gene"_hg38.bed \
  reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  expansion_report_"$gene" \
  alignments/????.bam
done
