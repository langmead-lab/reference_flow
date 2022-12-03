'''
The rules in this file performs coordinate system transformation for aligned reads
and sort the reads.
'''
rule lift_major_highq:
    input:
        sam = os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-mapqgeq{}.sam'.format(ALN_MAPQ_THRSD)),
        clft = os.path.join(DIR_MAJOR, EXP_LABEL + '-major.clft')
    output:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-mapqgeq{}-liftover.sam'.format(ALN_MAPQ_THRSD))
    params:
        os.path.join(DIR_FIRST_PASS, EXP_LABEL + '-major-mapqgeq{}-liftover'.format(ALN_MAPQ_THRSD))
    threads: THREADS
    run:
        shell('{LEVIOSAM2} lift -a {input.sam} -C {input.clft} -p {params} -t {threads}')

#: Refflow -- second pass
rule lift_refflow_secondpass_and_merge:
    input:
        maj_fp = os.path.join(
            DIR_FIRST_PASS,
            EXP_LABEL + '-major-mapqgeq{}-liftover.sam'.format(ALN_MAPQ_THRSD)),
        second_sam_path = os.path.join(
            DIR_SECOND_PASS,
            EXP_LABEL + '-major-{}-{}.merge_paths'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        clft_pop = expand(os.path.join(
            DIR_POP_GENOME, POP_DIRNAME + '/' +
            EXP_LABEL + '-' + POP_LEVEL + '_{GROUP}_' + POP_DIRNAME + '.clft'),
            GROUP = GROUP),
        clft_maj = os.path.join(DIR_MAJOR, EXP_LABEL + '-major.clft'),
    output:
        lfted_refflow_sam = os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        lfted_major_second_sam = os.path.join(DIR_SECOND_PASS, '2ndpass-maj-liftover.sam'),
        lfted_group_second_sam = [
            os.path.join(DIR_SECOND_PASS, '2ndpass-') + 
            g + '-liftover.sam' for g in GROUP]
    threads: THREADS
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

        # Copy lifted first pass sam.
        shell('cp {input.maj_fp} {output.lfted_refflow_sam};')
        # Run levioSAM and merge files.
        for i, sam in enumerate(list_sam):
            prefix = os.path.join(DIR,
                'experiments/' + wildcards.INDIV + '/' + POP_DIRNAME + 
                '/2ndpass-{}-liftover'.format(list_group[i]))
            if list_group[i] == 'maj':
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, input.lft_maj))
                shell('{LEVIOSAM2} lift -a {sam} -C {input.clft_maj} -p {prefix} -t {threads};')
                # Append reads to all-in-one lifted SAM.
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')
            elif list_group[i] in GROUP:
                for lft in input.clft_pop:
                    pop = os.path.basename(lft)
                    if lft.count(list_group[i]) > 0:
                        break
                sys.stderr.write('sam={}, lft = {}\n'.format(sam, lft))
                shell('{LEVIOSAM2} lift -a {sam} -C {clft} -p {prefix} -t {threads};')
                # Append reads to all-in-one lifted SAM.
                shell('grep -hv "^@" {prefix}.sam >> {output.lfted_refflow_sam};')

rule check_leviosam:
    input:
        expand(os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME)),
        INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'leviosam.done')))

'''
Sort SAM records
'''
rule sort_refflow:
    input:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover.sam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    output:
        os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted.bam'.format(ALN_MAPQ_THRSD, POP_DIRNAME))
    threads: 4
    run:
        shell('{SAMTOOLS} sort -@ {threads} -o {output} -O BAM {input};')

rule check_sort:
    input:
        expand(os.path.join(DIR_SECOND_PASS,
            EXP_LABEL + '-refflow-{}-{}-liftover-sorted.bam'.format(
                ALN_MAPQ_THRSD, POP_DIRNAME)), INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'sort.done')))

