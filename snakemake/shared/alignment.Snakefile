''' Reference flow: first pass '''
rule align_to_major:
    input:
        reads1 = READS1,
        reads2 = READS2,
        idx = expand(
            os.path.join(DIR, 'major/indexes/wg-maj.{idx}.bt2'),
            idx = IDX_ITEMS)
    params:
        index = os.path.join(DIR, 'major/indexes/wg-maj')
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major.sam')
    threads: THREADS
    shell:
        'bowtie2 --threads {THREADS} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}'
        #'bowtie2 --threads {THREADS} -x {params.index} -U {input.reads1} -S {output.sam}'

''' Refflow '''
rule refflow_split_aln_by_mapq:
    '''
    This is under a paired-end configuration
    '''
    input:
        sam = os.path.join(dir_first_pass, 'wg-major.sam')
    output:
        highq = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
        lowq = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD + '.sam'),
        lowq_reads1 = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
        lowq_reads2 = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_2.fq')
    params:
        fastq = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD),
        single_end = False,
        split_strategy = 'optimistic'
    shell:
        '{PYTHON} {DIR_SCRIPTS}/split_sam_by_mapq.py -s {input.sam} \
        -oh {output.highq} -ol {output.lowq} -oq {params.fastq} \
        -t {ALN_MAPQ_THRSD} --split-strategy {params.split_strategy}'

rule refflow_align_secondpass_paired_end:
    input:
        reads1 = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
        reads2 = os.path.join(DIR_FIRST_PASS,
            'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_2.fq'),
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
        'bowtie2 --reorder --threads {threads} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam};'

rule refflow_merge_secondpass:
    input:
        sam = expand(
            PREFIX_SECOND_PASS + '.sam',
            INDIV = INDIV, GROUP = GROUP),
        maj = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt{}-paired.sam'.format(ALN_MAPQ_THRSD))
    output:
        path = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        label = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.ids'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        merge_paths = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    params:
        prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
    run:
        # indiv = wildcards.INDIV
        dir_2p = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME)
        #maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/wg-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')

        #shell('echo {maj_lowq} > {output.path};')
        shell('echo {input.maj} > {output.path};')
        shell('echo "maj" > {output.label};')
        for g in GROUP:
            fn = os.path.join(dir_2p, 'wg-major-{}-{}-{}.sam'.format(ALN_MAPQ_THRSD, g, POP_DIRNAME))
            shell('ls {fn} >> {output.path};')
            shell('echo {g} >> {output.label};')
        shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
            -ids {output.label} -rs {RAND_SEED} -p {params.prefix} \
            -l {output.merge_paths} --paired-end')

'''
Paired-end, but a segment can align to different references, which is not accepted by GATK
'''
# rule refflow_split_aln_by_mapq:
#     input:
#         sam = os.path.join(dir_first_pass, 'wg-major.sam')
#     output:
#         highq = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqgeq' + ALN_MAPQ_THRSD + '.sam'),
#         lowq = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '.sam'),
#     shell:
#         'awk -v var="{ALN_MAPQ_THRSD}" \
#         \'{{ if ($5 >= var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
#         {output.highq};'
#         'awk -v var={ALN_MAPQ_THRSD} \
#         \'{{ if ($5 < var || $1 ~ /^@/) {{ print }} }}\' {input.sam} > \
#         {output.lowq};'
# 
# # split lowq sam into pairs and singletons
# rule refflow_split_lowq_sam:
#     input:
#         sam = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10.sam')
#     output:
#         singleton = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-singleton.sam'),
#         paired = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-paired.sam')
#     run:
#         fs_out = open(output.singleton, 'w')
#         f_out = open(output.paired, 'w')
#         with open(input.sam, 'r') as f:
#             name = ''
#             prev_line = ''
#             for line in f:
#                 if line[0] == '@':
#                     f_out.write(line)
#                     fs_out.write(line)
#                 else:
#                     new_name = line.split()[0]
#                     if new_name == name:
#                         name = ''
#                         f_out.write(prev_line)
#                         f_out.write(line)
#                         prev_line = ''
#                     elif name != '':
#                         fs_out.write(prev_line)
#                         name = new_name
#                         prev_line = line
#                     else:
#                         name = new_name
#                         prev_line = line
#             if prev_line:
#                 fs_out.write(prev_line)
# 
# rule refflow_sam_to_fastq_paired:
#     input:
#         paired = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-paired.sam')
#     output:
#         lowq_reads1 = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
#         lowq_reads2 = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_2.fq')
#     shell:
#         #'samtools fastq -1 {output.lowq_reads1} -2 {output.lowq_reads2} -n {output.lowq}'
#         'samtools fastq -n {input.paired} -1 {output.lowq_reads1} -2 {output.lowq_reads2}'
# 
# rule refflow_sam_to_fastq_singleton:
#     input:
#         singleton = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-singleton.sam'),
#     output:
#         lowq_reads_singleton = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_s.fq')
#     shell:
#         'samtools fastq -n {input.singleton} > {output.lowq_reads_singleton}'
# 
# rule refflow_align_secondpass_paired_end:
#     input:
#         reads1 = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_1.fq'),
#         reads2 = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_2.fq'),
#         #reads = os.path.join(DIR_FIRST_PASS,
#         #    'wg-major-mapqlt' + ALN_MAPQ_THRSD + '.fq'),
#         idx1 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.1.bt2'),
#         idx2 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.2.bt2'),
#         idx3 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.3.bt2'),
#         idx4 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.4.bt2'),
#         idx5 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.1.bt2'),
#         idx6 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.2.bt2')
#     params:
#         index = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX)
#     output:
#         sam = PREFIX_SECOND_PASS + '.sam'
#     threads: THREADS
#     shell:
#         'bowtie2 --reorder --threads {threads} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam};'
#         #'bowtie2 --reorder --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'
# 
# rule refflow_align_secondpass_singleton:
#     input:
#         reads = os.path.join(DIR_FIRST_PASS,
#             'wg-major-mapqlt' + ALN_MAPQ_THRSD + '_s.fq'),
#         idx1 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.1.bt2'),
#         idx2 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.2.bt2'),
#         idx3 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.3.bt2'),
#         idx4 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.4.bt2'),
#         idx5 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.1.bt2'),
#         idx6 = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX + '.rev.2.bt2')
#     params:
#         index = os.path.join(DIR_POP_GENOME_BLOCK_IDX, WG_POP_GENOME_SUFFIX)
#     output:
#         sam = PREFIX_SECOND_PASS + '-singleton.sam'
#     threads: THREADS
#     shell:
#         #'bowtie2 --reorder --threads {threads} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam};'
#         'bowtie2 --reorder --threads {threads} -x {params.index} -U {input.reads} -S {output.sam};'
# 
# rule refflow_merge_secondpass:
#     input:
#         sam = expand(
#             PREFIX_SECOND_PASS + '.sam',
#             INDIV = INDIV, GROUP = GROUP),
#         maj = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-paired.sam')
#     output:
#         path = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#         id = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.ids'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#         merge_paths = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
#     params:
#         prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
#     run:
#         # indiv = wildcards.INDIV
#         dir_2p = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME)
#         #maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/wg-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')
# 
#         #shell('echo {maj_lowq} > {output.path};')
#         shell('echo {input.maj} > {output.path};')
#         shell('echo "maj" > {output.id};')
#         for g in GROUP:
#             fn = os.path.join(dir_2p, 'wg-major-{}-{}-{}.sam'.format(ALN_MAPQ_THRSD, g, POP_DIRNAME))
#             shell('ls {fn} >> {output.path};')
#             shell('echo {g} >> {output.id};')
#         shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
#             -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
#             -l {output.merge_paths};')
# 
# rule refflow_merge_secondpass_singleton:
#     input:
#         sam = expand(
#             PREFIX_SECOND_PASS + '-singleton.sam',
#             INDIV = INDIV, GROUP = GROUP),
#         maj = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqlt10-singleton.sam')
#     output:
#         path = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}-singleton.paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#         id = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}-singleton.ids'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#         merge_paths = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}-singleton.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
#     params:
#         prefix = os.path.join(DIR_SECOND_PASS, '2ndpass')
#     run:
#         # indiv = wildcards.INDIV
#         dir_2p = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME)
#         #maj_lowq = os.path.join(DIR, 'experiments/' + wildcards.INDIV + '/wg-major-mapqlt' + ALN_MAPQ_THRSD + '.sam')
# 
#         #shell('echo {maj_lowq} > {output.path};')
#         shell('echo {input.maj} > {output.path};')
#         shell('echo "maj-singleton" > {output.id};')
#         for g in GROUP:
#             fn = os.path.join(dir_2p, 'wg-major-{}-{}-{}-singleton.sam'.format(ALN_MAPQ_THRSD, g, POP_DIRNAME))
#             shell('ls {fn} >> {output.path};')
#             shell('echo {g}-singleton >> {output.id};')
#         shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
#             -ids {output.id} -rs {RAND_SEED} -p {params.prefix} \
#             -l {output.merge_paths};')

rule check_alignment_refflow:
    input:
        merge_paths = expand(
            os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
            INDIV = INDIV),
#         merge_paths_singleton = expand(
#             os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}-singleton.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
#             INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'alignment_refflow.done')))
