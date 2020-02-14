'''
The rules in this file performs coordinate system transformation for aligned reads
and sort the reads.
'''

rule liftover_serialize_major:
    input:
        vcf_major = os.path.join(DIR, 'major/wg-maj.vcf'),
        length_map = LENGTH_MAP
    output:
        lft = os.path.join(DIR_MAJOR, 'wg-major.lft')
    params:
        os.path.join(DIR_MAJOR, 'wg-major')
    shell:
        '{LIFTOVER} serialize -v {input.vcf_major} -p {params} -k {input.length_map}'

rule liftover_serialize_pop_genome:
    input:
        vcf = os.path.join(DIR_POP_GENOME, POP_DIRNAME + '/' +
            'wg-superpop_{GROUP}_' + POP_DIRNAME  + '.vcf'),
        length_map = LENGTH_MAP
    output:
        lft = os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            'wg-superpop_{GROUP}_' + POP_DIRNAME + '.lft')
    params:
        os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            'wg-superpop_{GROUP}_' + POP_DIRNAME)
    run:
        shell('{LIFTOVER} serialize -v {input.vcf} -p {params} -k {input.length_map}')

rule liftover_lift_major_highq:
    input:
        sam = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqgeq{}.sam'.format(ALN_MAPQ_THRSD)),
        lft = os.path.join(DIR_MAJOR, 'wg-major.lft'),
        vcf_major = os.path.join(DIR, 'major/wg-maj.vcf')
    output:
        os.path.join(DIR_FIRST_PASS, 'wg-major-mapqgeq{}-liftover.sam'.format(ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, 'wg-major-mapqgeq{}-liftover'.format(ALN_MAPQ_THRSD))
    run:
        shell('{LIFTOVER} lift -a {input.sam} -l {input.lft} -p {params}')

#: Refflow -- second pass
rule liftover_lift_refflow_secondpass_and_merge:
    input:
        maj_fp = os.path.join(DIR_FIRST_PASS, 'wg-major-mapqgeq{}-liftover.sam'.format(ALN_MAPQ_THRSD)),
        second_sam_path = os.path.join(DIR_SECOND_PASS, 'wg-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        lft_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            'wg-superpop_{GROUP}_' + POP_DIRNAME + '.lft'),
            GROUP = GROUP),
        vcf_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            'wg-superpop_{GROUP}_' + POP_DIRNAME  + '.vcf'),
            GROUP = GROUP),
        lft_maj = os.path.join(DIR_MAJOR, 'wg-major.lft'),
        vcf_major = os.path.join(DIR_MAJOR, 'wg-maj.vcf')
    output:
        lfted_refflow_sam = os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        lfted_major_second_sam = os.path.join(DIR_SECOND_PASS, '2ndpass-maj-liftover.sam'),
        lfted_group_second_sam = [
            os.path.join(DIR_SECOND_PASS, '2ndpass-') + 
            g + '-liftover.sam' for g in GROUP]
    run:
        list_sam = []
        list_group = []
        #: files should be 
        #: DIR + '/experiments/{INDIV}/{POP_DIRNAME}/2ndpass-{}.sam'
        #: where g should be {GROUP} + 'maj'
        with open(input.second_sam_path, 'r') as f:
            for line in f:
                list_sam.append(line.rstrip())
                bn = os.path.basename(line)
                split_bn = os.path.splitext(bn)
                list_group.append(split_bn[0].split('-')[-1])
        for i, s in enumerate(list_sam):
            sys.stderr.write('sam={}, group = {}\n'.format(s, list_group[i]))
        
        #: copy lifted first pass sam 
        shell('cp {input.maj_fp} {output.lfted_refflow_sam};')
        for i in range(len(list_sam)):
            sam = list_sam[i]
            prefix = os.path.join(DIR,
                'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME + 
                '/2ndpass-{}-liftover'.format(list_group[i]))
            if list_group[i] == 'maj':
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, input.lft_maj))
                shell('{LIFTOVER} lift -a {sam} -l {input.lft_maj} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')
            elif list_group[i] in GROUP:
                for lft in input.lft_pop:
                    pop = os.path.basename(lft)
                    if lft.count(list_group[i]) > 0:
                        break
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, lft))
                shell('{LIFTOVER} lift -a {sam} -l {lft} -p {prefix};')
                #: append reads to all-in-one lifted SAM
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')

rule sort_refflow:
    input:
        os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    threads: 4
    run:
        shell('samtools sort -@ {threads} -o {output} -O BAM {input};')

rule check_sort:
    input:
        expand(os.path.join(DIR_SECOND_PASS,
            'wg-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)), INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'sorting.done')))

