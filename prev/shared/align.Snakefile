'''
Alignment using reference flow
    first pass reference genome: global major-allele reference
    second pass reference genome: population-specific references
'''

rule align_to_onekg_major:
    input:
        reads = PREFIX_PER + '_1.fq',
        idx = expand(PREFIX_MAJOR_IDX + '.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS)
    params:
        PREFIX_MAJOR_IDX
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major.sam'.format(CHROM))
    threads: THREADS
    shell:
        'bowtie2 --threads {threads} -x {params} -U {input.reads} -S {output.sam}'

rule refflow_separate_fistpass_results_onekg:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-major.sam'.format(CHROM))
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-mapqgeq{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        lowq_reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD))
    shell:
        'awk -v var="{ALN_MAPQ_THRSD}" \
        \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.highq};'
        'awk -v var={ALN_MAPQ_THRSD} \
        \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.lowq};'
        'samtools fastq {output.lowq} > {output.lowq_reads}'

rule refflow_align_secondpass_onekg:
    input:
        reads = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-mapqlt{}.fq'.format(CHROM, ALN_MAPQ_THRSD)),
        idx1 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.1.bt2',
        idx2 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.2.bt2',
        idx3 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.3.bt2',
        idx4 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.4.bt2',
        idx5 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.1.bt2',
        idx6 = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX + '.rev.2.bt2'
    params:
        index = DIR_POP_GENOME_BLOCK_IDX + POP_GENOME_SUFFIX
    output:
        sam = PREFIX_SECOND_PASS + '.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'

rule refflow_merge_secondpass_onekg:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '.sam',
            INDIV = INDIV, GROUP = GROUP),
        lowq = os.path.join(DIR_FIRST_PASS,
            'chr{}-major-mapqlt{}.sam'.format(CHROM, ALN_MAPQ_THRSD)),
    output:
        path = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        id = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}.ids'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
    run:
        indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + indiv + '/' + POP_DIRNAME)
        # maj_lowq = os.path.join(DIR, 'experiments/' + indiv + '/' + CHROM +
        #     '-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        shell('echo {input.lowq} > {output.path};')
        # shell('echo {maj_lowq} > {output.path};')
        shell('echo "major" > {output.id};')
        for g in GROUP:
            fn = os.path.join(dir_2p, 'chr{0}-major-{1}-{2}-{3}.sam'.format(CHROM, ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.id};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths};')

rule check_secondpass:
    input:
        # sam = expand(
        #     os.path.join(DIR_SECOND_PASS, '{0}-major-{1}-{2}.paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV, GROUP = GROUP
        # ),
        # gnomad = expand(
        #     os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}-gnomad.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
        #     INDIV = INDIV, GROUP = GROUP
        # ),
        onekg = expand(
            os.path.join(DIR_SECOND_PASS, 'chr{0}-major-{1}-{2}.merge_paths'.format(CHROM, ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV, GROUP = GROUP
        )
    output:
        touch(temp(os.path.join(DIR, 'refflow_secondpass.done')))
