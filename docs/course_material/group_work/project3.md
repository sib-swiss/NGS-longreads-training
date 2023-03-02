
## :material-bacteria: Project 3: Assembly and annotation of bacterial genomes

You will be working with PacBio sequencing data of five different bacterial strains. Divide the strains over the members of the group and generate an assembly and annotation.

!!! info "Project aim"
    Generate and evaluate an assembly of a bacterial genome out of PacBio reads. 

There are five different strains: 



Each strain has a tarfile available. Download only the data for the strains that you will require: 

```sh
mkdir -p ~/workdir/groupwork_assembly
cd ~/workdir/groupwork_assembly

# change this to your strain:
STRAIN="LWX12"

wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/group_work_assembly/"$STRAIN".tar.gz
tar -xvf "$STRAIN".tar.gz
rm "$STRAIN".tar.gz
```

The downloaded directory has the following structure (here's an example for LWH7):

```

```

### Before you start

You can start this project with dividing the strains over the different group members. In principle, each group member will go through all the steps of assembly and annotation:

1. Quality control with `NanoPlot`
2. Assembly with `flye`
3. Assembly QC with `BUSCO`
4. Annotation with `prokka`

You can do this for both the CLR reads and HiFi reads and compare the results. 

### Tasks and questions

!!! note
    You have **four** cores available. Use them! For most tools you can specificy the number of cores/cpus as an argument. 

!!! note
    All require software can be found in the conda environment `assembly`. Load it like this:

    ```sh
    conda activate assembly
    ```

    **Before you run `prokka`**

    The `conda` installation misses a perl module. Install it in the `assembly` environment like this:

    ```sh
    cpanm Bio::SearchIO::hmmer --force
    ```

* Perform a quality control with `NanoPlot`.
    * How is the read quality? Is this quality expected?
    * How is the read length?
* Perform an assembly with `flye`. 
    * Have a look at the helper first with `flye --help`. Make sure you pick the correct mode (i.e. `--pacbio-??`). 
    * Check out the output. Where is the assembly? How is the quality? For that, check out `assembly_info.txt`. 
* Check the completeness with `BUSCO`. Have a good look at the manual first. You can use automated lineage selecton by specifying `--auto-lineage-prok`. After you have run `BUSCO`, you can generate a nice completeness plot with `generate_plot.py`. You can check its usage with `generate_plot.py --help`. 
    * How is the completeness? Is this expected?
* Perform an annotation with `prokka`. Again, check the manual first. After the run, have a look at for example the statistics in `PROKKA_[date].txt`. For a nice table of annotated genes have a look in `PROKKA_[data].tsv`. 



* Compare the assembly and annotation between the Illumina, CLR and HiFi reads. Do you see any differences? 
* Compare the assemblies of the different strains. Are assembly qualities similar? Can you think of reasons why?
* **BONUS**: Polish the CLR assembly with the Illumina reads by using `pilon`. For this you will need to align the Illumina reads to the assembly first. Use `minimap2` for that while setting `-x` to `sr`. For pilon, specify the resulting bam file by using the option `--frags`. 
    * Does the polishing improve the assembly? Why (not)?
