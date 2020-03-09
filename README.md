# Reference Flow

Reference flow is a computational genomics pipeline to align short sequencing 
reads with a number of reference genomes built with population genomic information included.

The preprint ["Reducing reference bias using multiple population reference genomes"](https://doi.org/10.1101/2020.03.03.975219) is available on bioRxiv.

The experiments we've done are provided in the [Reference Flow Experiments repository](https://github.com/langmead-lab/reference_flow-experiments).

## Preparation

### Install Snakemake

Reference flow is written using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and it can be installed using conda:

```
conda install -c bioconda -c conda-forge snakemake
```

Other installation approaches are also provided in the [Snakemake installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Download the 1000 Genomes Project population table

[Population information for the 1000 Genomes Project individuals](https://www.internationalgenome.org/faq/which-samples-are-you-sequencing/) can be downloaded by running the script:

```
sh src/download_1kg_pop_table.sh
```

### Liftover

Liftover is used to perform coordinate system mappign between population reference genomes and a standard reference genome such as GRCh38. 
The liftover software tracks the "edits" between a pair of genomes using a VCF file.
Please refer to the [liftover github page](https://github.com/alshai/liftover) to install the software.


## Running reference flow with pre-built indexes

```
sh src/download_prebuilt_indexes.sh
cd snakemake
snakemake -j 32
```

Before executing the Snakemake pipeline, it is recommended to perform a dry-run `snakemake -np` to preview the scheduling plan. The option `-j` specifies the number of threads used, which is set to 32 in the above example.

By default, a directory called `run` will be created under `snakemake` and all the results will be under it. 
The alignment results are aggrelated as a single SAM file, which uses the GRCh38 coordinate system by default.
The final SAM file will appear in directory `snakemake/run/experiments/test/thrds0_S1_b1000_ld1/wg-refflow-10-thrds0_S1_b1000_ld1-liftover.sam`.

To use another read set, please change the `READS1` parameter in `snakemake/config.yaml`, 
or run command:

```
snakemake -j 32 --config READS1=<file>
```

## Running complete reference flow pipeline

By set the `USE_PREBUILT` option to `False`, users can run the complete reference flow pipeline.

### Download Reference genome

The GRCh38 reference genome can be downloaded from NCBI using the following script:

```
sh src/download_genome.sh
```

Users may choose other reference genomes of interest by changing the `GENOME` parameter in `snakemake/config.yaml`.

### Download GRCh38 call sets from the 1000 Genomes Project

[The 1000 Genomes Project GRCh38 call sets with SNVs and indels](https://www.internationalgenome.org/announcements/Variant-calls-from-1000-Genomes-Project-data-on-the-GRCh38-reference-assemlby/) can be downloaded using script:

```
sh src/download_1kg_vcf.sh
```

Users may switch to the 1KG phase-3 call set by using the corresponding VCFs and changing parameters `DIR_VCF`, `VCF_PREFIX`, and `VCF_SUFFIX` in the configuration file.
For now we only tested on call sets provided by the 1000 Genomes Project. 
Call sets provided by other studies may differ slightly in population labelling and VCF format, causing the pipeline not working properly.

## Setting reference flow configuration

Parameters for reference flow are specified in `snakemake/config.yaml` and described [here](snakemake/README.md).

We recommend users try a dry-run by `snakemake -np` to check if the pipeline is ready.

When all set, run

```
snakemake -j 32
```
