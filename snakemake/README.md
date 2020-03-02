# Reference Flow parameters

- `READS1` : reads to align

- `INDIV` : name of the tested sample. This parameter is used to name directories and files.

- `DIR` : directory where the outputs will be put

- `THREADS` : max number of threads used for each rule.
This is different from the number specified by `snakemake -j <t>`, which is the number of threads for the entire Snakemake pipeline. I.e., there can be multiple programs running with `THREADS` until the sum reaches `<t>`.

- `GENOME` : reference genome; usually a standard GRC genome

- `ALN_MAPQ_THRSD` : mapping quality cutoff to split read into committed and deferred groups

- `RAND_SEED`: random seed used in the reference flow stochastic reference genome update process and for aligner

### Population dataset parameters

- `DIR_VCF` : directory where the VCFs are put

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

