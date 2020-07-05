rule build_major:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf.gz'),
        genome = GENOME
    output:
        vcfgz = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz'),
        vcfgz_idx = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz.csi'),
        out_genome = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa'),
    shell:
        '{BCFTOOLS} view -O z -q 0.5000001 -G -o {output.vcfgz} -v snps,indels -m2 -M2 {input.vcf};'
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
