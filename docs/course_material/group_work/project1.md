
## :fontawesome-solid-brain: Project 1: Differential isoform expression analysis of ONT data

In this project, you will be working with data from the same resource as the data we have already worked on:

> Padilla, Juan-Carlos A., Seda Barutcu, Ludovic Malet, Gabrielle Deschamps-Francoeur, Virginie Calderon, Eunjeong Kwon, and Eric Lécuyer. “Profiling the Polyadenylated Transcriptome of Extracellular Vesicles with Long-Read Nanopore Sequencing.” BMC Genomics 24, no. 1 (September 22, 2023): 564. https://doi.org/10.1186/s12864-023-09552-6.

It is Oxford Nanopore Technology sequencing data of cDNA from extracellular vesicles and whole cells. It is primarily used to discover new splice variants. We will use the dataset to do that and in addition do a differential isoform expression analysis with [FLAIR](https://github.com/BrooksLabUCSC/flair).

!!! info "Project aim"
    Discover new splice variants and identify differentially expressed isoforms.

You can download the required data like this:

```sh
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz
```

!!! note
    Download the data file package in your shared working directory, i.e. : `/group_work/<group name>`. Only one group member has to do this. You can add the group work directory to the workspace in VScode by opening the menu on the top right (hamburger symbol), click **File** > **Add folder to workspace** and type the path to the group work directory.

This will create a directory `project1` with the following structure:

```
project1/
├── reads
│   ├── Cell_1.fastq.gz
│   ├── Cell_2.fastq.gz
│   ├── Cell_3.fastq.gz
│   ├── EV_1.fastq.gz
│   ├── EV_2.fastq.gz
│   └── EV_3.fastq.gz
├── reads_manifest.tsv
└── references
    ├── Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf
    └── Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa

2 directories, 9 files
```

In the reads folder a fastq file with reads, which are described in `reads_manifest.csv`. EV means 'extracellular vesicle', Cell means 'entire cells'. In the references folder you can find the reference sequence and annotation.

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

* Have a look at the [FLAIR documentation](https://flair.readthedocs.io/en/latest/index.html).
* FLAIR and all its dependencies are in the the pre-installed conda environment named `flair`. You can activate it with `conda activate flair`.
* Merge the separate alignments with `samtools merge`, index the merged bam file, and generate a `bed12` file with the command `bam2Bed12`
* Run `flair correct` on the `bed12` file. Add the `gtf` to the options to improve the alignments.
* Run `flair collapse` to generate isoforms from corrected reads. This steps takes ~1.5 hours to run.
* Generate a count matrix with `flair quantify` by using the isoforms fasta and `reads_manifest.tsv` (takes ~45 mins to run).

!!! danger "Paths in `reads_manifest.tsv`"
    The paths in `reads_manifest.tsv` are relative, e.g. `reads/striatum-5238-batch2.fastq.gz` points to a file relative to the directory from which you are running `flair quantify`. So the directory from which you are running the command should contain the directory `reads`. If not, modify the paths in the file accordingly (use full paths if you are not sure).

* Now you can do several things:
    * Do a differential expression analysis. In `scripts/` there's a basic R script to do the analysis. Go to your specified IP and port to login to RStudio server (the username is `rstudio`).
    * Investigate the isoform usage with the flair script `plot_isoform_usage.py`
    * Investigate productivity of the different isoforms.