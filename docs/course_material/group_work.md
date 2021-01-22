
The last part of this course will consist of project-based-learning. This means that you will work in groups on a single question. We will split up into groups of five people.

!!! note "If working independently"
    If you are working independently, you probably can not work in a group. However, you can test your skills with these real biological datasets. Realize that the datasets and calculations are (much) bigger compared to the exercises, so check if your computer is up for it. You'll probably need around 4 cores, 16G of RAM and 10G of harddisk.

!!! note "If online"
    If the course takes place online, we will use break-out rooms to communicate within groups. Please stay in the break-out room during the day, also if you are working individually.

## Roles & organisation

Project based learning is about learning by doing, but also about *peer instruction*. This means that you will be both a learner and a teacher. There will be differences in levels among participants, but because of that, some will learn efficiently from people that have just learned, and others will teach and increase their understanding.

Each project has **tasks** and **questions**. By performing the tasks, you should be able to answer the questions. At the start of the project, make sure that each of you gets a task assigned. You should consider the tasks and questions as a guidance. If interesting questions pop up during the project, you are **encouraged** to work on those. Also, you don't have to perform all the tasks and answer all the questions.

In the afternoon of day 1, you will divide the initial tasks, and start on the project. On day 2, you can work on the project in the morning and in the first part of the afternoon. We will conclude the projects with a **10-minute presentation** of each group.

## :fontawesome-solid-brain: Project 1: Differential isoform expression analysis of ONT data

In this project, you will be working with data from the same resource as the data we have already worked on:

Clark, M. B. et al (2020). *Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain*. Molecular Psychiatry, 25(1), 37–47. [https://doi.org/10.1038/s41380-019-0583-1](https://doi.org/10.1038/s41380-019-0583-1).

It is Oxford Nanopore Technology sequencing data of PCR amplicons of the gene CACNA1C. It is primarily used to discover new splice variants. We will use the dataset to do that and in addition do a differential isoform expression analysis with [FLAIR](https://github.com/BrooksLabUCSC/flair).

You can download the required data like this:

```sh
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/groupwork_ont.tar.gz
tar -xvf groupwork_ont.tar.gz
rm groupwork_ont.tar.gz
```

This will create a directory `groupwork_ont` with the following structure:

```
groupwork_ont
├── alignments
│   ├── cerebellum-5238-batch2.bam
│   ├── cerebellum-5298-batch2.bam
│   ├── cerebellum-5346-batch2.bam
│   ├── parietal_cortex-5238-batch1.bam
│   ├── parietal_cortex-5298-batch1.bam
│   └── parietal_cortex-5346-batch1.bam
├── counts
│   └── counts_matrix_test.tsv
├── reads
│   ├── cerebellum-5238-batch2.fastq.gz
│   ├── cerebellum-5298-batch2.fastq.gz
│   ├── cerebellum-5346-batch2.fastq.gz
│   ├── parietal_cortex-5238-batch1.fastq.gz
│   ├── parietal_cortex-5298-batch1.fastq.gz
│   ├── parietal_cortex-5346-batch1.fastq.gz
│   ├── striatum-5238-batch2.fastq.gz
│   ├── striatum-5298-batch2.fastq.gz
│   └── striatum-5346-batch2.fastq.gz
└── scripts
    └── differential_expression_example.Rmd

4 directories, 17 files
```

Download the fasta file and gtf like this:

```sh
cd groupwork_ont/
mkdir reference
wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
wget ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
gunzip *.gz
```

### Before you start

You can start this project with dividing initial tasks. Because some intermediate files are already given, participants can develop scripts/analyses at different steps of the full analysis from the start. Possible starting points are:

* Quality control, running `fastqc` and `NanoPlot`
* Alignment, running `minimap2`
* Develop scripts required to run FLAIR
* Differential expression analysis.

### Tasks & questions

* Perform QC with `fastqc` and with `NanoPlot`.
    * Do you see a difference between them?
    * How is the read quality compared to the publication?
* Align each sample separately with `minimap2` with default parameters. Set parameters `-x` and `-G` to the values we have used during the [QC and alignment exercises](../qc_alignment#3-read-alignment). You can use 4 threads (set the number of threads with `-t`)

!!! danger "Start the alignment on day 1"
    The alignment takes about 6 minutes per sample, so in total about one hour to run. Make sure you have started the alignment at day 1, so you don't have to wait for the results on day 2. Use `nohup myscript.sh &` to be able to logout while `myscript.sh` is running (`tmux` and `screen` are also available).

* Clone the [FLAIR repository](https://github.com/BrooksLabUCSC/flair) to the server, and check out the documentation.
* Merge the separate alignments with `samtools merge`, index the merged bam file, and generate a `bed12` file with the script `flair/bin/bam2Bed12.py`
* Run `flair.py correct` on the `bed12` file. Add the `gtf` to the options to improve the alignments.
* Run `flair.py collapse` to generate isoforms from corrected reads. This steps takes ~1 hour to run.
* Generate a count matrix with `flair.py quantify` by using the isoforms fasta and `reads_manifest.tsv`.

!!! danger "Paths in `reads_manifest.tsv`"
    The paths in `reads_manifest.tsv` are relative, e.g. `reads/striatum-5238-batch2.fastq.gz` points to a file relative to the directory from which you are running `flair.py quantify`. So the directory from which you are running the command should contain the directory `reads`. If not, modify the paths in the file accordingly (use full paths if you are not sure).

* Now you can do several things:
    * Do a differential expression analysis. In `scripts/` there's a basic R script to do the analysis. Go to `[SERVERIP]:8787` to login to RStudio server.
    * Investigate the isoform usage with the flair script `plot_isoform_usage.py`
    * Investigate productivity of the different isoforms.

## :fontawesome-solid-disease: Project 2: Repeat expansion analysis of PacBio data

You will be working with data from an experiment in which DNA of 8 individuals was sequenced for five different targets by using Pacbio's no-Amp targeted sequencing system. Two of these targets contain repeat expansions that are related to a disease phenotype.

| individual 	| disease1 	| disease2 	|
|------------	|----------	|----------	|
| 1015       	| disease  	| healthy  	|
| 1016       	| disease  	| healthy  	|
| 1017       	| disease  	| healthy  	|
| 1018       	| disease  	| healthy  	|
| 1019       	| healthy  	| healthy  	|
| 1020       	| healthy  	| disease  	|
| 1021       	| healthy  	| disease  	|
| 1022       	| healthy  	| disease  	|

You can get the reads and sequence targets with:

```sh
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/groupwork_pacbio.tar.gz
tar -xvf groupwork_pacbio.tar.gz
rm groupwork_pacbio.tar.gz
```

It has the following directory structure:

```
groupwork_pacbio
├── alignments
│   ├── bc1020.aln.bam
│   ├── bc1021.aln.bam
│   └── bc1022.aln.bam
├── reads
│   ├── 1015.fastq.gz
│   ├── 1016.fastq.gz
│   ├── 1017.fastq.gz
│   ├── 1018.fastq.gz
│   ├── 1019.fastq.gz
│   ├── 1020.fastq.gz
│   ├── 1021.fastq.gz
│   └── 1022.fastq.gz
└── targets
    ├── target_gene1_hg38.bed
    └── target_gene2_hg38.bed

3 directories, 13 files
```

The targets in gene1 and gene2 are described in `targets/target_gene1_hg38.bed` and `targets/target_gene2_hg38.bed` respectively. The columns in these `.bed` files describe the chromosome, start, end, name, motifs, and whether the motifs are in reverse complement.

You can download the reference genome like this:

```sh
cd groupwork_pacbio
mkdir reference
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Before you start

You can start this project with dividing initial tasks. Because some intermediate files are already given, participants can develop scripts/analyses at different steps of the full analysis from the start. Possible starting points are:

* Browse IGV to find the genes
* Perform the QC with `NanoPlot`
* Perform the alignment with `minimap2`
* Do the repeat analysis with `makeReports.sh`

Alignment files to do an initial repeat analysis are in the `tar.gz` package. However, it contains only the files for individuals with disease2. You can develop scripts and analyses based on that. To do the full analysis, all the alignments will need to be run.

### Tasks & questions

* Load the bed files into IGV and navigate to the regions they annotate.
    * In which genes are the targets?
    * What kind of diseases are associated with these genes?
* Perform a quality control with `NanoPlot`.
    * How is the read quality? These are circular concensus sequences (ccs). Is this quality expected?
    * How is the read length?
* Align the reads to hg38 with `minimap2`. For the option `-x` you can use `asm20`. Generate separate alignment files for each individual.

!!! note "Alternatively use [`pbmm2`](https://github.com/PacificBiosciences/pbmm2)"
    Pacific Biosciences has developed a wrapper for `minimap2` that contains settings specific for PacBio reads, named [`pbmm2`](https://github.com/PacificBiosciences/pbmm2). It might slightly improve your alignments. It is installed in the conda environment. Feel free to give it a try if you have time left.

* Clone the PacBio [apps-scripts repository](https://github.com/PacificBiosciences/apps-scripts.git) to the server. The script apps-scripts/RepeatAnalysisTools/makeReports.sh generates repeat expansion reports. Check out the [documentation](https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools#auto-generate-all-reports), and generate repeat expansion reports for all individuals on both gene1 and gene2.
* Check out the report output and read the further [documentation of RepeatAnalysisTools](https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools#auto-generate-all-reports).
    - How is the enrichment?
    - Does the clustering make sense? How does the clustering look in IGV?
    - Which individual is affected with which disease?
    - Based on the size of the expansions, can you say something about expected disease severity?

> This tutorial is based on data provided by Pacific Biosciences at https://downloads.pacbcloud.com/public/dataset/RepeatExpansionDisorders_NoAmp/
