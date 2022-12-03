rule build_major:
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf.gz'),
        genome = GENOME
    output:
        vcf = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf'),
        vcfgz = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz'),
        vcfgz_idx = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.vcf.gz.csi'),
        out_genome = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa'),
        out_genome_fai = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa.fai'),
        inverted_chain = temp(os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.inverted.chain'))
    shell:
        '{BCFTOOLS} view -O z -q 0.5000001 -G -o {output.vcf} -v snps,indels -m2 -M2 {input.vcf};'
        '{BGZIP} -c {output.vcf} > {output.vcfgz};'
        '{BCFTOOLS} index {output.vcfgz};'
        '{BCFTOOLS} consensus -c {output.inverted_chain} -f {input.genome} -o {output.out_genome} {output.vcfgz};'
        '{SAMTOOLS} faidx {output.out_genome}'

rule invert_major_chain:
    input:
        inverted_chain = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.inverted.chain')
    output:
        chain = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.chain')
    shell:
        '{PYTHON} {CHAINTOOLS}/src/invert.py -c {input.inverted_chain} -o {output.chain}'

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

