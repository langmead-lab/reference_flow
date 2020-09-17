''' This snakefile includes rules that perform single-end alignment.'''
''' Perform first-pass alignment.'''
rule align_to_major:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR, 'major/indexes/' + EXP_LABEL + '-maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'major/indexes/' + EXP_LABEL + '-maj')
    output:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

''' Split first-pass alignment into high-quality and low-quality files. '''
rule refflow_split_aln_by_mapq:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major.sam')
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-major-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        lowq = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-major-mapqlt' + ALN_MAPQ_THRSD + '.sam'),
        lowq_reads = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
    shell:
        'awk -v var="{ALN_MAPQ_THRSD}" \
        \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.highq};'
        'awk -v var={ALN_MAPQ_THRSD} \
        \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
        {output.lowq};'
        'samtools fastq {output.lowq} > {output.lowq_reads}'

''' Align low-quality reads using population genomes.'''
rule refflow_align_secondpass_single_end:
    input:
        reads1 = os.path.join(DIR_FIRST_PASS,
            EXP_LABEL + '-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
        idx1 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.1.bt2'),
        idx2 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.2.bt2'),
        idx3 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.3.bt2'),
        idx4 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.4.bt2'),
        idx5 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.1.bt2'),
        idx6 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.2.bt2')
    params:
        index = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX)
    output:
        sam = PREFIX_SECOND_PASS + '.sam'
    threads: THREADS
    shell:
        'bowtie2 --reorder --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

''' Merge refflow results. '''
rule refflow_merge_secondpass:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '.sam',
            INDIV = INDIV, GROUP = GROUP),
        maj = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-mapqlt{}.sam'.format(ALN_MAPQ_THRSD))
    output:
        path = os.path.join(
            DIR_SECOND_PASS,
            EXP_LABEL + '-major-{}-{}.paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        label = os.path.join(
            DIR_SECOND_PASS,
            EXP_LABEL + '-major-{}-{}.ids'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(
            DIR_SECOND_PASS,
            EXP_LABEL + '-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
    run:
        dir_2p = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME)
        shell('echo {input.maj} > {output.path};')
        shell('echo "maj" > {output.label};')
        for g in GROUP:
            fn = os.path.join(
                dir_2p,
                EXP_LABEL + '-major-{}-{}-{}.sam'.format(
                    ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.label};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.label} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths}')

rule check_alignment_refflow:
    input:
        merge_paths = expand(
            os.path.join(DIR_SECOND_PASS,
                         EXP_LABEL + '-major-{}-{}.merge_paths'.format(
                            ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_refflow.done')))
