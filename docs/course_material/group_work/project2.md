

## :fontawesome-solid-disease: Project 2: Repeat expansion analysis of PacBio data

You will be working with data from an experiment in which DNA of 8 individuals was sequenced for five different targets by using Pacbio's no-Amp targeted sequencing system. Two of these targets contain repeat expansions that are related to a disease phenotype.

!!! info "Project aim"
    Estimate variation in repeat expansions in two target regions, and relate them to a disease phenotype.

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

!!! note
    Download the data file package in your shared working directory, i.e. : `/group_work/<group name>` or `~/<group name>`. Only one group member has to do this.

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
cd reference
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

!!! danger "Start the alignment as soon as possible"
    The alignment takes quite some time. Try to start the alignment as soon as possible. You can speed up your alignment by first making an index, e.g.:

    ```sh
    minimap2 \
    -x asm20 \
    -d reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.mmi \
    reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```

    Refer to the generated index (`.mmi` file) as reference in the alignment command, e.g.:

    ```
    minimap2 \
    -a \
    -x asm20 \
    -t 4 \
    reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.mmi \
    reads/<my_reads.fastq.gz>
    ```

!!! note "Alternatively use [`pbmm2`](https://github.com/PacificBiosciences/pbmm2)"
    Pacific Biosciences has developed a wrapper for `minimap2` that contains settings specific for PacBio reads, named [`pbmm2`](https://github.com/PacificBiosciences/pbmm2). It might slightly improve your alignments. It is installed in the conda environment. Feel free to give it a try if you have time left.

* Clone the PacBio [apps-scripts repository](https://github.com/PacificBiosciences/apps-scripts.git) to the server. All dependencies are in the conda environment `pacbio`. Activate it with `conda activate pacbio`. The script apps-scripts/RepeatAnalysisTools/makeReports.sh generates repeat expansion reports. Check out the [documentation](https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools#auto-generate-all-reports), and generate repeat expansion reports for all individuals on both gene1 and gene2.
* Check out the report output and read the further [documentation of RepeatAnalysisTools](https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools#auto-generate-all-reports).
    - How is the enrichment?
    - Does the clustering make sense? How does the clustering look in IGV?
    - Which individual is affected with which disease?
    - Based on the size of the expansions, can you say something about expected disease severity?

> This tutorial is based on data provided by Pacific Biosciences at https://downloads.pacbcloud.com/public/dataset/RepeatExpansionDisorders_NoAmp/
