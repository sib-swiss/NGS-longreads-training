
cd ~/project/project1/

mkdir -p flair_align

flair align \
-g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa \
-r reads/*.fastq.gz \
--output flair_align/flair.aligned

mkdir -p flair_correct

flair correct -q flair_align/flair.aligned.bed \
-f references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf \
-g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa \
--output flair_correct/flair

mkdir -p flair_collapse

flair collapse \
-g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa \
-q flair_correct/flair_all_corrected.bed \
-r reads/*.fastq.gz \
--gtf references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf \
--output flair_collapse/flair.collapse

mkdir -p flair_quantify

flair quantify \
-r reads_manifest.tsv \
-i flair_collapse/flair.collapse.isoforms.fa \
--output flair_quantify/flair.quantify \
--sample_id_only

plot_isoform_usage \
flair_collapse/flair.collapse.isoforms.bed \
flair_quantify/flair.quantify.counts.tsv ENSG00000113013
