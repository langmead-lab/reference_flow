'''
Rules for building population genomes
'''
rule prepare_pop_indiv:
    output:
        expand(
            os.path.join(DIR, '1KG_indivs/sample_superpop_{GROUP}.txt'),
            GROUP = GROUP
        )
    params:
        prefix = os.path.join(DIR, '1KG_indivs/sample')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/list_indiv_from_pop.py '
        '-p {FAMILY} -sp {SPOP} -op {params.prefix}'

rule build_pop_vcf:
    input:
        vcf = PREFIX_VCF_F + '.vcf',
        indiv_group = os.path.join(
            DIR,
            '1KG_indivs/sample_superpop_{GROUP}.txt'
        )
    output:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}.vcf.gz')
    shell:
        '{BCFTOOLS} view -S {input.indiv_group} '
        '--force-samples {input.vcf} -V mnps,other -m2 -M2 | '
        'bgzip > {output.vcf_gz}'

rule get_pop_sample:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}.vcf.gz')
    output:
        vcf_header = os.path.join(DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}.samples')
    shell:
        '{BCFTOOLS} view -h {input.vcf_gz} | tail -1 '
        '> {output.vcf_header}'

rule filter_pop_vcf:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}.vcf.gz'),
        vcf_header = os.path.join(DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}.samples')
    output:
        vcf = os.path.join(
            DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}_t' + str(POP_THRSD) + '.vcf'
        )
    run:
        fn = list({input.vcf_header})[0]
        with open(fn, 'r') as f:
            for line in f:
                n = len(line.split()) - 9
                thrsd = int(n * 2 * float(POP_THRSD))
                filt = 'AC > {}'.format(thrsd)
                break
        shell('{BCFTOOLS} view -i "{filt}" \
            -v snps,indels {input.vcf_gz} > {output.vcf};')

rule build_pop_genome:
    input:
        vcf = os.path.join(
            DIR_POP_GENOME,
            'chr{CHROM}_superpop_{GROUP}_t' + str(POP_THRSD) + '.vcf'
        )
    output:
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.fa'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.var'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.vcf'
        )
    params:
        prefix = os.path.join(
            DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX
        )
    run:
        if POP_STOCHASTIC == 1 and POP_USE_LD == 1:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --chrom {wildcards.CHROM} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --stochastic -rs {RAND_SEED} \
                --block-size {POP_BLOCK_SIZE} --ld')
        elif POP_STOCHASTIC == 1:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --chrom {wildcards.CHROM} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --stochastic -rs {RAND_SEED} \
                --block-size {POP_BLOCK_SIZE}')
        else:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --chrom {wildcards.CHROM} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels')

rule merge_pop_genome:
    input:
        expand(os.path.join(DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.fa'),
            CHROM = CHROM, GROUP = GROUP)
    output:
        os.path.join(DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.fa')
    run:
        list_fn = []
        for fn in input:
            print (fn)
            print (wildcards.GROUP)
            if fn.count(wildcards.GROUP) > 0:
                list_fn.append(fn)
        shell('cat {list_fn} >> {output}')

rule aggregate_pop_vcf:
    input:
        vcf = expand(os.path.join(DIR_POP_GENOME_BLOCK,
            POP_GENOME_SUFFIX + '.vcf'),
            CHROM = CHROM, GROUP = GROUP)
    output:
        vcf = os.path.join(DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.vcf')
    run:
        list_fn = []
        for fn in input:
            print (fn)
            print (wildcards.GROUP)
            if fn.count(wildcards.GROUP) > 0:
                list_fn.append(fn)
        shell('{BCFTOOLS} concat -o {output.vcf} {list_fn}')

rule build_pop_genome_index:
    input:
        genome = os.path.join(DIR_POP_GENOME_BLOCK, WG_POP_GENOME_SUFFIX + '.fa')
    output:
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.1.bt2'),
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.2.bt2'),
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.3.bt2'),
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.4.bt2'),
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.1.bt2'),
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.2.bt2')
    params:
        os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX)
    threads: THREADS
    shell:
        'bowtie2-build --threads {threads} {input.genome} {params};'

rule check_pop_genome:
    input:
        expand(
            DIR_POP_GENOME_BLOCK_IDX + WG_POP_GENOME_SUFFIX + '.{IDX_ITEMS}.bt2',
            GROUP = GROUP, IDX_ITEMS = IDX_ITEMS
        )
    output:
        touch(temp(os.path.join(DIR, 'prepare_pop_genome.done')))
