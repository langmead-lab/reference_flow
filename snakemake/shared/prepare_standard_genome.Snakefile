# rule build_major:
#     input:
#         vcf = PREFIX_VCF_F + '.vcf',
#         genome_chrom = os.path.join(DIR, 'chr{CHROM}.fa')
#     output:
#         vcf_major = PREFIX_MAJOR_F + '.vcf',
#         out_genome = PREFIX_MAJOR + '.fa',
#         out_var = PREFIX_MAJOR + '.var',
#         out_vcf = PREFIX_MAJOR + '.vcf'
#     params:
#         out_prefix = os.path.join(DIR, 'major/chr{CHROM}_maj')
#     shell:
#         # '{BCFTOOLS} view -O z -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
#         # '{output.vcf_major_gz};'
#         # '{BCFTOOLS} index {output.vcf_major_gz};'
#         # 'bgzip -cd {output.vcf_major_gz} > {output.vcf_major};'
#         '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
#         '{output.vcf_major};'
#         '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
#         '    --ref {input.genome_chrom} --vcf {output.vcf_major} '
#         '    --chrom {wildcards.CHROM} --out-prefix {params.out_prefix} '
#         '    --include-indels'
# 
# rule merge_major_fasta:
#     input:
#         expand(PREFIX_MAJOR + '.fa', CHROM = CHROM)
#     output:
#         os.path.join(DIR, 'major/wg-maj.fa')
#     shell:
#         'cat {input} >> {output}'
# 
# rule aggregate_major_vcf:
#     input:
#         vcf = expand(PREFIX_MAJOR + '.vcf', CHROM = CHROM)
#     output:
#         vcf = os.path.join(DIR, 'major/wg-maj.vcf')
#     shell:
#         '{BCFTOOLS} concat -o {output.vcf} {input.vcf}'

# rule build_major:
#     input:
#         vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf')
#     output:
#         vcf_major = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj-bcftools.vcf'),
#         out_genome = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa'),
#         out_var = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.var'),
#         out_vcf = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf')
#     params:
#         out_prefix = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj')
#     shell:
#         '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
#         '{output.vcf_major};'
#         '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
#         '   --ref {GENOME} --vcf {output.vcf_major} '
#         '   --out-prefix {params.out_prefix} '
#         '   --include-indels'

rule build_major:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf'),
        genome = GENOME
    output:
        vcf = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf'),
        vcfgz = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz'),
        vcfgz_idx = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz.csi'),
        out_genome = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa'),
    shell:
        '{BCFTOOLS} view -O v -q 0.5000001 -G -o {output.vcf} -v snps,indels -m2 -M2 {input.vcf};'
        'bgzip -c {output.vcf} > {output.vcfgz};'
        '{BCFTOOLS} index {output.vcfgz};'
        '{BCFTOOLS} consensus -f {input.genome} -o {output.out_genome} {output.vcfgz}'

rule build_major_index:
    input:
        os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa')
    output:
        expand(
            os.path.join(DIR, 'major/indexes/' + EXP_LABEL + '-maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'major/indexes/' + EXP_LABEL + '-maj')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'
 
rule check_standard_genomes:
    input:
        expand(
            os.path.join(DIR, 'major/indexes/' + EXP_LABEL + '-maj.{idx}.bt2'),
            idx = IDX_ITEMS),
    output:
        touch(temp(os.path.join(DIR, 'prepare_standard_genome.done')))
