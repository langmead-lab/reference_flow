'''
Rules for building population genomes
'''
rule prepare_pop_indiv:
    output:
        expand(
            os.path.join(DIR, '1KG_indivs/sample_' + POP_LEVEL + '_{GROUP}.txt'), GROUP = GROUP)
    params:
        prefix = os.path.join(DIR, '1KG_indivs/sample')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/list_indiv_from_pop.py '
        '-p {FAMILY} -sp {SPOP} -op {params.prefix}'

rule build_pop_vcf:
    '''
    Filter VCF by population groups.
    Each output VCF includes only indivs in the specified population group.
    '''
    input:
        vcf = os.path.join(DIR, EXP_LABEL + '_filtered.vcf.gz'),
        indiv_group = os.path.join(DIR, '1KG_indivs/sample_' + POP_LEVEL + '_{GROUP}.txt')
    output:
        vcf_gz = os.path.join(DIR_POP_GENOME, EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}.vcf.gz')
    shell:
        '{BCFTOOLS} view -S {input.indiv_group} '
        '--force-samples {input.vcf} -V mnps,other -m2 -M2 | bgzip > {output.vcf_gz}'

rule get_pop_sample:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}.vcf.gz')
    output:
        vcf_header = os.path.join(DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}.samples')
    shell:
        '{BCFTOOLS} view -h {input.vcf_gz} | tail -1 '
        '> {output.vcf_header}'

rule filter_pop_vcf:
    input:
        vcf_gz = os.path.join(DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}.vcf.gz'),
        vcf_header = os.path.join(DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}.samples')
    output:
        vcf = os.path.join(
            DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}_t' + str(POP_THRSD) + '.vcf'
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

rule construct_pop_vcf:
    input:
        vcf = os.path.join(
            DIR_POP_GENOME,
            EXP_LABEL + '_' + POP_LEVEL + '_{GROUP}_t' + str(POP_THRSD) + '.vcf'
        )
    output:
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.var'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.vcf.gz'
        ),
        os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.vcf.gz.tbi'
        )
    params:
        prefix = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX
        )
    run:
        if POP_STOCHASTIC == 1 and POP_USE_LD == 1:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --stochastic -rs {RAND_SEED} \
                --block-size {POP_BLOCK_SIZE} --ld --var-only')
        elif POP_STOCHASTIC == 1:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --stochastic -rs {RAND_SEED} \
                --block-size {POP_BLOCK_SIZE} --var-only')
        else:
            shell('{PYTHON} {DIR_SCRIPTS}/update_genome.py \
                --ref {GENOME} --vcf {input.vcf} \
                --out-prefix {params.prefix} \
                --include-indels --var-only')
        shell('{BGZIP} {params.prefix}.vcf')
        shell('{TABIX} {params.prefix}.vcf.gz')

rule build_pop_genome:
    input:
        vcf = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.vcf.gz'
        )
    output:
        fasta = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.fa'
        ),
        fai = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.fa.fai'
        ),
        inverted_chain = temp(os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.inverted.chain'
        )),
    shell:
        '{BCFTOOLS} consensus -f {GENOME} -c {output.inverted_chain} {input.vcf} > {output.fasta};'
        '{SAMTOOLS} faidx {output.fasta}'

rule invert_pop_chain:
    input:
        inverted_chain = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.inverted.chain'
        )
    output:
        chain = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.chain'
        )
    shell:
        '{PYTHON} {CHAINTOOLS}/src/invert.py -c {input.inverted_chain} -o {output.chain}'

rule leviosam2_index_chain_major:
     input:
         chain = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.chain'),
         fai = os.path.join(DIR, 'major/' + EXP_LABEL + '-maj.fa.fai'),
     output:
         clft = os.path.join(DIR_MAJOR, EXP_LABEL + '-major.clft')
     params:
         prefix = os.path.join(DIR_MAJOR, EXP_LABEL + '-major')
     shell:
         '{LEVIOSAM2} index -c {input.chain} -F {input.fai} -p {params.prefix}'

rule leviosam2_index_chain_pop:
    input:
        chain = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.chain'
        ),
        fai = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.fa.fai'
        ),
    output:
        clft = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX + '.clft'
        )
    params:
        prefix = os.path.join(
            DIR_POP_GENOME_BLOCK,
            WG_POP_GENOME_SUFFIX
        )
    shell:
        '{LEVIOSAM2} index -c {input.chain} -F {input.fai} -p {params.prefix};'

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
            GROUP=GROUP, IDX_ITEMS=IDX_ITEMS, POP_LEVEL=POP_LEVEL
        ),
    output:
        touch(temp(os.path.join(DIR, 'prepare_pop_genome.done')))

