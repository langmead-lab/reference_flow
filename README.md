# Reference flow

Reference flow is a computational genomics pipeline to align short sequencing 
reads with a number of reference genomes built with population genomic information included.

The preprint will be posted on bioRxiv.


## Preparation

### Insall Snakemake

Reference flow is written using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and it can be installed using conda:

```
conda install -c bioconda -c conda-forge snakemake
```

Other installation approaches are also provided in the [Snakemake installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Reference genome

The GRCh38 reference genome can be downloaded from NCBI using the following script. Users may choose other reference genomes of interest.

```
sh src/download_genome.sh
```

### GRCh38 call sets from the 1000 Genomes Project

#todo

### Download 1000 Genomes Project population table

[Population information for 1000 Genomes Project individuals](https://www.internationalgenome.org/faq/which-samples-are-you-sequencing/) can be downloaded by running

```
sh resources/download_1kg_pop_table.sh
```

### Liftover

#todo

## Reference flow configuration

Parameters for reference flow are specified in `snakemake/config.yaml`.

## Running reference flow

Example:

```
cd snakemake
snakemake -j 32
```

By default, a directory called `run` will be created under `snakemake` and all the results will be under it. Final alignment outputs will appear in directory `run/experiments/<indiv>/<refflow_config>/`.
