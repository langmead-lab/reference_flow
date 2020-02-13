# Reference flow

Reference flow is a computational genomics pipeline to align short sequencing 
reads with a number of reference genomes built with population genomic information included.

The preprint will be posted on bioRxiv.


## Preparation

### Install Snakemake

Reference flow is written using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and it can be installed using conda:

```
conda install -c bioconda -c conda-forge snakemake
```

Other installation approaches are also provided in the [Snakemake installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Download Reference genome

The GRCh38 reference genome can be downloaded from NCBI using the following script:

```
sh src/download_genome.sh
```

Users may choose other reference genomes of interest by changing the `GENOME` parameter in `snakemake/config.yaml`.

### GRCh38 call sets from the 1000 Genomes Project

#todo

### Download 1000 Genomes Project population table

[Population information for 1000 Genomes Project individuals](https://www.internationalgenome.org/faq/which-samples-are-you-sequencing/) can be downloaded by running the script:

```
sh resources/download_1kg_pop_table.sh
```

### Liftover

Liftover is used to perform coordinate system mappign between population reference genomes and a standard reference genome such as GRCh38. 
The liftover software tracks the "edits" between a pair of genomes using a VCF file.
Please refer to the [liftover github page](https://github.com/alshai/liftover) to install the software.

## Reference flow configuration

Parameters for reference flow are specified in `snakemake/config.yaml`.

## Running reference flow

Example:

```
cd snakemake
snakemake -j 32
```

Before executing the Snakemake pipeline, it is recommended to perform a dry-run `snakemake -np` to preview the scheduling plan. The option `-j` specifies the number of threads used, which is set to 32 in the above example.

By default, a directory called `run` will be created under `snakemake` and all the results will be under it. 
The alignment results are aggrelated as a single SAM file, which is sorted and all using the GRCh38 coordinate system by default.
The final SAM file will appear in directory `run/experiments/<indiv>/<refflow_config>/`.

