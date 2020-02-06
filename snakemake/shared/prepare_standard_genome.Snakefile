rule build_major:
    input:
        vcf = PREFIX_VCF_F + '.vcf',
        genome_chrom = os.path.join(DIR, 'chr{CHROM}.fa')
    output:
        vcf_major = PREFIX_MAJOR_F + '.vcf',
        out_genome = PREFIX_MAJOR + '.fa',
        out_var = PREFIX_MAJOR + '.var',
        out_vcf = PREFIX_MAJOR + '.vcf'
    params:
        out_prefix = os.path.join(DIR, 'major/chr{CHROM}_maj')
    shell:
        # '{BCFTOOLS} view -O z -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
        # '{output.vcf_major_gz};'
        # '{BCFTOOLS} index {output.vcf_major_gz};'
        # 'bgzip -cd {output.vcf_major_gz} > {output.vcf_major};'
        '{BCFTOOLS} view -O v -q 0.5 {input.vcf} -e \'AF = 0.5\' -v snps,indels -m2 -M2 > '
        '{output.vcf_major};'
        '{PYTHON} {DIR_SCRIPTS}/update_genome.py '
        '    --ref {input.genome_chrom} --vcf {output.vcf_major} '
        '    --chrom {wildcards.CHROM} --out-prefix {params.out_prefix} '
        '    --include-indels'

rule merge_major_fasta:
    input:
        expand(PREFIX_MAJOR + '.fa', CHROM = CHROM)
    output:
        os.path.join(DIR, 'major/wg-maj.fa')
    shell:
        'cat {input} >> {output}'

rule aggregate_major_vcf:
    input:
        vcf = expand(PREFIX_MAJOR + '.vcf', CHROM = CHROM)
    output:
        vcf = os.path.join(DIR, 'major/wg-maj.vcf')
    shell:
        '{BCFTOOLS} concat -o {output.vcf} {input.vcf}'

rule build_major_index:
    input:
        os.path.join(DIR, 'major/wg-maj.fa')
    output:
        expand(
            os.path.join(DIR, 'major/indexes/wg-maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        os.path.join(DIR, 'major/indexes/wg-maj')
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input} {params}'
 
rule check_standard_genomes:
    input:
        expand(
            os.path.join(DIR, 'major/indexes/wg-maj.{idx}.bt2'),
            idx = IDX_ITEMS),
    output:
        touch(temp(os.path.join(DIR, 'prepare_standard_genome.done')))
