
## :material-bacteria: Project 3: Assembly and annotation of bacterial genomes

You will be working with PacBio sequencing data of eight different bacterial species. Divide the species over the members of the group and generate an assembly and annotation. After that, guess the species. 

!!! info "Project aim"
    Generate and evaluate an assembly of a bacterial genome out of PacBio reads. 

There are eight different species: `sample_[1-8].fastq.gz` 



Each species has a fastq file available. You can download all fastq files like this: 

```sh
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project3.tar.gz
tar -xvf project3.tar.gz
rm project3.tar.gz
```

!!! note
    Download the data file package in your shared working directory, i.e. : `/group_work/<group name>` or `~/<group name>`. Only one group member has to do this.

This will create a directory `project3` with the following structure:

```
project3
|-- sample_1.fastq.gz
|-- sample_2.fastq.gz
|-- sample_3.fastq.gz
|-- sample_4.fastq.gz
|-- sample_5.fastq.gz
|-- sample_6.fastq.gz
|-- sample_7.fastq.gz
`-- sample_8.fastq.gz

0 directories, 8 files
```

### Before you start

You can start this project with dividing the species over the different group members. In principle, each group member will go through all the steps of assembly and annotation:

1. Quality control with `NanoPlot`
2. Assembly with `flye`
3. Assembly QC with `BUSCO`
4. Annotation with `prokka`

### Tasks and questions

!!! note
    You have **four** cores available. Use them! For most tools you can specificy the number of cores/cpus as an argument. 

!!! note
    All require software can be found in the conda environment `assembly`. Load it like this:

    ```sh
    conda activate assembly
    ```

* Perform a quality control with `NanoPlot`.
    * How is the read quality? Is this quality expected?
    * How is the read length?
* Perform an assembly with `flye`. 
    * Have a look at the helper first with `flye --help`. Make sure you pick the correct mode (i.e. `--pacbio-??`). 
    * Check out the output. Where is the assembly? How is the quality? For that, check out `assembly_info.txt`. 
    * What species did you assemble? Choose from this list:
    ```
    Acinetobacter baumannii
    Bacillus cereus
    Bacillus subtilis
    Burkholderia cepacia
    Burkholderia multivorans
    Enterococcus faecalis
    Escherichia coli
    Helicobacter pylori
    Klebsiella pneumoniae
    Listeria monocytogenes
    Methanocorpusculum labreanum
    Neisseria meningitidis
    Rhodopseudomonas palustris
    Salmonella enterica
    Staphylococcus aureus
    Streptococcus pyogenes
    Thermanaerovibrio acidaminovorans
    Treponema denticola
    Vibrio parahaemolyticus
    ```
    * Did flye assemble any plasmid sequences?
* Check the completeness with `BUSCO`. Have a good look at the manual first. You can use automated lineage selecton by specifying `--auto-lineage-prok`. After you have run `BUSCO`, you can generate a nice completeness plot with `generate_plot.py`. You can check its usage with `generate_plot.py --help`. 
    * How is the completeness? Is this expected?
* Perform an annotation with `prokka`. Again, check the manual first. After the run, have a look at for example the statistics in `PROKKA_[date].txt`. For a nice table of annotated genes have a look in `PROKKA_[date].tsv`. 
* Compare the assemblies of the different species. Are assembly qualities similar? Can you think of reasons why?

> This tutorial is based on data provided by Pacific Biosciences at https://downloads.pacbcloud.com/public/dataset/2021-11-Microbial-96plex/