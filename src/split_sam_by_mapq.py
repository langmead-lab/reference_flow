import argparse
import threading

class ReadWriteLock:
    ''' From https://www.oreilly.com/library/view/python-cookbook/0596001673/ch06s04.html '''
    """ A lock object that allows many simultaneous "read locks", but
    only one "write lock." """

    def __init__(self):
        self._read_ready = threading.Condition(threading.Lock())
        self._readers = 0

    def acquire_read(self):
        """ Acquire a read lock. Blocks only if a thread has
        acquired the write lock. """
        self._read_ready.acquire()
        try:
            self._readers += 1
            #print ('acquire_read()', self._readers)
        finally:
            self._read_ready.release()
            #print ('acquire_read(): _read_ready.release()')

    def release_read(self):
        """ Release a read lock. """
        self._read_ready.acquire()
        try:
            self._readers -= 1
            if not self._readers:
                self._read_ready.notifyAll()
        finally:
            self._read_ready.release()

    def acquire_write(self):
        """ Acquire a write lock. Blocks until there are no
        acquired read or write locks. """
        self._read_ready.acquire()
        while self._readers > 0:
            self._read_ready.wait()

    def release_write(self):
        """ Release a write lock. """
        self._read_ready.release()


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
        '-pe', '--paired-end', action='store_true',
        help="Set if reads are paired-end [Off]"
    )
    parser.add_argument(
        '-ss', '--split-strategy',
        help='Split strategy (only needed when using paired-end reads): \
                "optimistic/opt" takes the pair if any of it is high quality, \
                "pessimistic/pes" takes the pair when both are high quality'
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

def get_reverse_complement(seq):
    '''
    Perform reverse and complement transformation
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def write_line_to_fastq(line, f_fastq):
    '''
    Write a SAM line to a FASTQ file
    '''
    fields = line.split()
    is_reverse = int(fields[1]) & 16
    if is_reverse:
        # reverse and complement SEQ, reverse QUAL
        f_fastq.write('@' + fields[0] + '\n')
        f_fastq.write(get_reverse_complement(fields[9]) + '\n')
        f_fastq.write('+\n')
        f_fastq.write(fields[10][::-1] + '\n')
    else:
        f_fastq.write('@' + fields[0] + '\n')
        f_fastq.write(fields[9] + '\n')
        f_fastq.write('+\n')
        f_fastq.write(fields[10] + '\n')

def process_paired_end_data_line(
    line,
    line_nxt,
    fhigh_out,
    flow_out,
    flow_fq1_out,
    flow_fq2_out,
    mapq_threshold,
    split_strategy
):
    fields = line.split()
    fields_nxt = line_nxt.split()
    try:
        assert(fields[0] == fields_nxt[0])
    except:
        print ('Warning: singleton read')
        print (line.rstrip())
        print (line_nxt.rstrip())
        exit(1)

    if split_strategy in ['opt', 'optimistic']:
        mapq = max(int(fields[4]), int(fields_nxt[4]))
    else:
        mapq = min(int(fields[4]), int(fields_nxt[4]))

    concordant = (fields[2] == fields_nxt[2])

    if mapq >= mapq_threshold and concordant:
        fhigh_out.write(line)
        fhigh_out.write(line_nxt)
    else:
        flow_out.write(line)
        flow_out.write(line_nxt)
        if flow_fq1_out:
            flag = int(fields[1])
            flag_nxt = int(fields_nxt[1])
            # first/second segment: current/next line
            if (flag & 64) and (flag_nxt & 128):
                write_line_to_fastq(line, flow_fq1_out)
                write_line_to_fastq(line_nxt, flow_fq2_out)
            # first/second segment: next/current line
            elif (flag & 128) and (flag_nxt & 64):
                write_line_to_fastq(line, flow_fq2_out)
                write_line_to_fastq(line_nxt, flow_fq1_out)
            else:
                print ('Error: read is not paired-end')
                print (line)
                print (line_nxt)
                exit (1)

def process_paired_end_data_parallel_core(
    f_in, fhigh_out, flow_out, flow_fq1_out, flow_fq2_out, mapq_threshold, split_strategy, rw_lock
):  
    chunk_size = 10000
    lines = []
    while 1:
        # spins until the lock is released
        rw_lock.acquire_write()

        # reads header. Headers are not counted in the chunk.
        while 1:
            line = f_in.readline()
            if line and (line[0] == '@'):
                fhigh_out.write(line)
                flow_out.write(line)
            else:
                # non-header line or end-of-file
                break
        if line:
            lines.append(line)
        for i in range(chunk_size - 1):
            line_nxt = f_in.readline()
            if line_nxt:
                lines.append(line_nxt)
        rw_lock.release_write()

        if not lines:
            return

        for i in range(0, len(lines), 2):
            process_paired_end_data_line(lines[i], lines[i+1], fhigh_out, flow_out, flow_fq1_out, flow_fq2_out, mapq_threshold, split_strategy)
        lines = []

def process_paired_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold, split_strategy):
    '''
    Process paired-end data
    '''
    if fastq_prefix:
        flow_fq1_out = open(fastq_prefix + '_1.fq', 'w')
        flow_fq2_out = open(fastq_prefix + '_2.fq', 'w')

    # rw_lock = ReadWriteLock()
    # threads = []
    # for i in range(4):
    #     t = threading.Thread(
    #         target = process_paired_end_data_parallel_core,
    #         args=(f_in, fhigh_out, flow_out, flow_fq1_out, flow_fq2_out, mapq_threshold, split_strategy, rw_lock)
    #     )
    #     t.start()
    #     threads.append(t)

    # for t in threads:
    #     print (t)
    #     t.join()

    for line in f_in:
        if line[0] == '@':
            fhigh_out.write(line)
            flow_out.write(line)
        else:
            line_nxt = f_in.readline()
            process_paired_end_data_line(line, line_nxt, fhigh_out, flow_out, flow_fq1_out, flow_fq2_out, mapq_threshold, split_strategy)

def split_sam_by_mapq(args):
    mapq_threshold = args.mapq_threshold
    is_paired_end = args.paired_end
    if not is_paired_end:
        split_strategy = ''
    else:
        split_strategy = args.split_strategy
        assert split_strategy in ['optimistic', 'opt', 'pessimistic', 'pes']

    f_in = open(args.sam, 'r')
    fhigh_out = open(args.sam_highq, 'w')
    flow_out = open(args.sam_lowq, 'w')
    fastq_prefix = args.fastq_lowq_prefix
    
    if not is_paired_end:
        process_single_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold)
    else:
        process_paired_end_data(f_in, fhigh_out, flow_out, fastq_prefix, mapq_threshold, split_strategy)
    
if __name__ == '__main__':
    args = parse_args()
    split_sam_by_mapq(args)
