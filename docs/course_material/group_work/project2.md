

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
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project2.tar.gz
tar -xvf project2.tar.gz
rm project2.tar.gz
```

!!! note
    Download the data file package in your shared working directory, i.e. : `/group_work/<group name>`. Only one group member has to do this.

It has the following directory structure:

```
project2
├── reads
│   ├── 1015.fastq.gz
│   ├── 1016.fastq.gz
│   ├── 1017.fastq.gz
│   ├── 1018.fastq.gz
│   ├── 1019.fastq.gz
│   ├── 1020.fastq.gz
│   ├── 1021.fastq.gz
│   └── 1022.fastq.gz
├── reference
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa
│   └── Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa.fai
└── targets
    └── targets.bed

3 directories, 11 files
```

The targets in gene1 and gene2 are described in `targets/targets.bed`. The columns in these `.bed` files describe the chromosome, start, end, and describe the motifs. To reduce computational load, the reference contains only chromosome 4 and X of the hg38 human reference genome. 

### Tasks & questions

!!! info "Activate the conda environment"
    The tools you will be needed for these exercises are in the conda environment `lr-tools`. Every time you open a new terminal, activate it with:

    ```sh
    conda activate lr-tools
    ```

* Load the bed files into IGV and navigate to the regions they annotate.
    * In which genes are the targets?
    * What kind of diseases are associated with these genes?
* Perform a quality control with `NanoPlot`.
    * How is the read quality? These are circular concensus sequences (ccs). Is this quality expected?
    * How is the read length?
* Align the reads to `reference/Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa` with `minimap2`. For the option `-x` you can use `asm20`. Generate separate alignment files for each individual. Check out some of the bam files in IGV. How does that look? 

!!! note "Alternatively use [`pbmm2`](https://github.com/PacificBiosciences/pbmm2)"
    Pacific Biosciences has developed a wrapper for `minimap2` that contains settings specific for PacBio reads, named [`pbmm2`](https://github.com/PacificBiosciences/pbmm2). It might slightly improve your alignments. It is installed in the conda environment. Feel free to give it a try if you have time left.

* Use `trgt` to genotype the repeats. Basically, you want to know the expansion size of each repeat in each sample. Based on this, you can figure out which sample has abnormal expansions in which repeat. To run `trgt` read [the manual](https://github.com/PacificBiosciences/trgt/blob/main/docs/tutorial.md). After the alignment, all required input files should be there. 
* To visualize the output, use samtools to sort and index the bam file with the reads spanning the repeats (this is also explained in the manual - no need to run `bcftools`).
*  Run `trvz` to visualize the output. The allele plot should suffice. The visualization will give you a nice overview of the repeat expansions in the samples. Based on the different sizes of the repeat expansions, can you relate the repeat expansions to the disease phenotype?

> This tutorial is based on data provided by Pacific Biosciences at https://downloads.pacbcloud.com/public/dataset/RepeatExpansionDisorders_NoAmp/
