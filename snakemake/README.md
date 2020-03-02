# Reference Flow parameters

### Data-dependent parameters
For users interested in the default RandFlow-LD pipeline based on GRCh38 and 1KG GRCh38 call set,
changing parameters in this section is sufficient.

- `READS1` : reads to align

- `INDIV` : name of the tested sample. This parameter is used to name directories and files.

- `DIR` : directory where the outputs will be put

- `THREADS` : max number of threads used for each rule.
This is different from the number specified by `snakemake -j <t>`, which is the number of threads for the entire Snakemake pipeline. I.e., there can be multiple programs running with `THREADS` until the sum reaches `<t>`.


### Pipeline options
- `USE_PREBUILT` : whether to use the RandFlow-LD pre-built indexes

- `SORT_SAM` : whether to sort the final SAM output

- `ALN_MAPQ_THRSD` : mapping quality cutoff to split read into committed and deferred groups


### Dataset parameters
- `GENOME` : reference genome; usually a standard GRC genome

- `DIR_VCF` : directory where the VCFs are put

- `VCF_PREFIX`, `VCF_SUFFIX` : prefix and suffix of a VCF file.
For example, we set `VCF_PREFIX` to "ALL.chr" and 
`VCF_SUFFIX` to ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
for the 1KG GRCh38 call set, where the VCFs are named as 
"ALL.chr<chrom>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

- `CHROM` : a list specifying chromosomes included. 
Reference flow finds VCFs corresponding to the specified chromosomes and builds population genomes for them.

- `GROUP` : population groups included to build second pass population genomes


### Second-pass genome parameters
- `POP_THRSD` : allele frequency threshold. 0: do not filter by frequency; 0.5: only use major alleles
  
- `POP_STOCHASTIC` : whether to use stochastic update. 1: stochastic update; 0: deterministic

- `POP_BLOCK_SIZE` : size of phase-preserving LD blocks. 
Set to 1 when doing independent sampling. Only applicable to call sets providing phase information.

- `POP_USE_LD` : whether to include local LD information. 1: phase-preserving; 0: independent-sampling

We provided configurations for `MajorFlow`, `RandFlow` and `RandFlow-LD` in the file.
Users may uncomment the setting of interest.


### Individual-population mappings
- `FAMILY` : mapping between individuals and populations
- `SPOP` : mapping between populations and super populations

Currently we only support call sets from the 1000 Genomes Project, including phase-3 and GRCh38 call sets.
Default mappings are provided under `../resources` directory.
For users interested in using other call sets, 
changes in the pipeline may be needed depending on the labelling and VCF format adopted by the call set.
`FAMILY` and `SPOP` are used under the `prepare_pop_indiv` rule in `shared/prepare_pop_genome.Snakefile` and 
the `build_dict_pop_to_spop` rule in `shared/functions.Snakefile`. 


### Chromosomal mappings
- `LENGTH_MAP` : lengths for each chromosome
- `CHROM_MAP` : mapping from `1` to `chr1`, etc


### Paths of software
- `BCFTOOLS` : we used bcftools 1.9-206-g4694164. 
There is a `bcftools consensus` issue with the major release so we used a development version.

- `SAMTOOLS` : we used samtools v1.9

- `LIFTOVER` : `../liftover/liftover`

- `PYTHON`: we developed Reference flow using Python3.7

- `DIR_SCRIPTS` : `../src`

### Randomness control
- `RAND_SEED`: random seed used in the reference flow stochastic reference genome update process and for aligner
