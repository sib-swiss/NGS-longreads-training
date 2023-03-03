
## :fontawesome-solid-brain: Project 1: Differential isoform expression analysis of ONT data

In this project, you will be working with data from the same resource as the data we have already worked on:

Clark, M. B. et al (2020). *Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain*. Molecular Psychiatry, 25(1), 37–47. [https://doi.org/10.1038/s41380-019-0583-1](https://doi.org/10.1038/s41380-019-0583-1).

It is Oxford Nanopore Technology sequencing data of PCR amplicons of the gene CACNA1C. It is primarily used to discover new splice variants. We will use the dataset to do that and in addition do a differential isoform expression analysis with [FLAIR](https://github.com/BrooksLabUCSC/flair).

!!! info "Project aim"
    Discover new splice variants and identify differentially expressed isoforms.

You can download the required data like this:

```sh
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/groupwork_ont.tar.gz
tar -xvf groupwork_ont.tar.gz
rm groupwork_ont.tar.gz
```

!!! note
    Download the data file package in your shared working directory, i.e. : `/group_work/<group name>` or `~/<group name>`. Only one group member has to do this.

This will create a directory `groupwork_ont` with the following structure:

```
groupwork_ont/
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
├── reads_manifest.tsv
└── scripts
    └── differential_expression_example.Rmd

4 directories, 18 files
```

Download the fasta file and gtf like this:

```sh
cd groupwork_ont/
mkdir reference
cd reference
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

!!! danger "Start the alignment as soon as possible"
    The alignment takes about 6 minutes per sample, so in total about one hour to run. Try to start the alignment as soon as possible. You can speed up your alignment by first making an index, e.g.:

    ```sh
    minimap2 \
    -x splice \
    -d reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
    reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa
    ```

    Refer to the generated index (`.mmi` file) as reference in the alignment command, e.g.:

    ```
    minimap2 \
    -a \
    -x splice \
    -G 500k \
    -t 4 \
    reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
    reads/<my_reads.fastq.gz>
    ```

* Clone the [FLAIR repository](https://github.com/BrooksLabUCSC/flair) to the server, and check out the documentation. All FLAIR dependencies are in the the pre-installed conda environment named `flair`. You can activate it with `conda activate flair`.
* Merge the separate alignments with `samtools merge`, index the merged bam file, and generate a `bed12` file with the script `flair/bin/bam2Bed12.py`
* Run `flair correct` on the `bed12` file. Add the `gtf` to the options to improve the alignments.
* Run `flair collapse` to generate isoforms from corrected reads. This steps takes ~1 hour to run.
* Generate a count matrix with `flair quantify` by using the isoforms fasta and `reads_manifest.tsv`.

!!! danger "Paths in `reads_manifest.tsv`"
    The paths in `reads_manifest.tsv` are relative, e.g. `reads/striatum-5238-batch2.fastq.gz` points to a file relative to the directory from which you are running `flair quantify`. So the directory from which you are running the command should contain the directory `reads`. If not, modify the paths in the file accordingly (use full paths if you are not sure).

* Now you can do several things:
    * Do a differential expression analysis. In `scripts/` there's a basic R script to do the analysis. Go to your specified IP and port to login to RStudio server (the username is `rstudio`).
    * Investigate the isoform usage with the flair script `plot_isoform_usage.py`
    * Investigate productivity of the different isoforms.