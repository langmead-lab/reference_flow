import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--sam',
        help='Original SAM file'
    )
    parser.add_argument(
        '-oh', '--sam-highq',
        help='High quality alignments (in SAM format)'
    )
    parser.add_argument(
        '-ol', '--sam-lowq',
        help='Low quality alignments (in SAM format)'
    )
    parser.add_argument(
        '-oq', '--fastq-lowq-prefix',
        help='Predix of low quality reads (in FASTQ format). \
                There will be two files for a SAM with paired-end data. \
                Leave empty to not write to FASTQ [None]'
    )
    parser.add_argument(
        '-t', '--mapq-threshold', type=int, default=10,
        help='Mapping quality threshold, reads with MAPQ greater or \
                equal to `-t` are called high quality [10]'
    )
    parser.add_argument(
        '--single-end', action='store_true',
        help="Set if reads are single-end [Off]"
    )
    parser.add_argument(
        '--split-strategy',
        help='Split strategy (only needed when using paired-end reads): \
                "optimistic" takes the pair if any of it is high quality, \
                "pessimistic" takes the pair when both are high quality'
    )
    # parser.add_argument(
    #     '-rs', '--rand-seed',
    #     help='random seed for controlled randomness [None]'
    # )
    # parser.add_argument(
    #     '-p', '--prefix',
    #     help='prefix of the merged files'
    # )
    args = parser.parse_args()
    return args

def write_line_to_fastq(line, f_fastq):
    '''
    Write a SAM line to a FASTQ file
    '''
    fields = line.split()
    f_fastq.write('@' + fields[0] + '\n')
    f_fastq.write(fields[9] + '\n')
    f_fastq.write('+\n')
    f_fastq.write(fields[10] + '\n')

def process_single_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold):
    '''
    Process single-end data, don't need to worry about the strategies for paired-end reads
    '''
    f_fastq = fastq_prefix + '.fq'
    for line in f_in:
        if line[0] == '@':
            fhigh_out.write(line)
            flow_out.write(line)
        else:
            mapq = line.split()[4]
            if mapq >= mapq_threshold:
                fhigh_out.write(line)
            else:
                flow_out.write(line)
                write_line_to_fastq(line, f_fastq)

def process_paired_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold, split_strategy):
    '''
    Process paired-end data
    '''
    flow_fq1_out = open(fastq_prefix + '_1.fq', 'w')
    flow_fq2_out = open(fastq_prefix + '_2.fq', 'w')

    name = ''
    prev_line = ''
    prev_mapq = 0
    for line in f_in:
        if line[0] == '@':
            fhigh_out.write(line)
            flow_out.write(line)
        else:
            new_name = line.split()[0]
            if new_name == name:
                # see a pair
                if split_strategy == 'optimistic':
                    mapq = max(prev_line.split()[4], line.split()[4])
                else:
                    mapq = min(prev_line.split()[4], line.split()[4])

                if mapq >= mapq_threshold:
                    fhigh_out.write(prev_line)
                    fhigh_out.write(line)
                else:
                    flow_out.write(prev_line)
                    flow_out.write(line)
                    write_line_to_fastq(line, flow_fq1_out)
                    write_line_to_fastq(line, flow_fq2_out)
                name = ''
                prev_line = ''
                prev_mapq = 0
            elif (new_name != name) and (name == ''):
                name = new_name
                prev_line = line
                prev_mapq = line[4]
            else:
                # singleton, should not happen
                print ('Error: read {} is a singleton. Please check the data'.format(name))
                exit(1)
    if prev_line:
        print ('Error: read {} is a singleton. Please check the data'.format(name))
        exit(1)


def split_sam_by_mapq(args):
    mapq_threshold = args.mapq_threshold
    is_single_end = args.single_end
    if is_single_end:
        split_strategy = ''
    else:
        split_strategy = args.split_strategy
        assert split_strategy in ['optimistic', 'pessimistic']

    f_in = open(args.sam, 'r')
    fhigh_out = open(args.sam_highq, 'w')
    flow_out = open(args.sam_lowq, 'w')
    fastq_prefix = args.fastq_lowq_prefix
    
    if is_single_end:
        process_single_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold)
    else:
        process_paired_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold, split_strategy)
    
if __name__ == '__main__':
    args = parse_args()
    split_sam_by_mapq(args)
