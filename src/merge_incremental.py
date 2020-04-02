import sys
import random
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-ns', '--sam-list',
        help='list of paths of SAM files'
    )
    parser.add_argument(
        '-ids', '--id-list',
        help='list of ids of files'
    )
    parser.add_argument(
        '-l', '--log',
        help='log file to write paths of merged files'
    )
    parser.add_argument(
        '-rs', '--rand-seed',
        help='random seed for controlled randomness [None]'
    )
    parser.add_argument(
        '-p', '--prefix',
        help='prefix of the merged files'
    )
    parser.add_argument(
        '--paired-end', action='store_true',
        help="Set if reads are paired-end [Off]"
    )
    args = parser.parse_args()
    return args

def get_info_from_sam_line(line, flag = False):
    '''
    Read one line from SAM and return some info.
    Default: QNAME, alignment score (AS:i) and MAPQ
    Optional: flag

    Args:
        line: a line from SAM
        flag: set True to return flag [False]

    Returns:
        QNAME (string)
        AS:i (int)
        MAPQ (int)
        flag (int) if `flag == True`
    '''
    if line[0] == '@':
        return 'header', None
    line = line.split()
    name = line[0]
    flag = int(line[1])
    mapq = int(line[4])
    score = 1 #: represents unmapped
    for i in line:
        if i.startswith('AS:i'):
            score = int(i.split(':')[-1])

    if not flag:
        return name, score, mapq
    else:
        return name, score, mapq, flag

def compare_score_and_mapq(list_info):
    '''
    Select the best record given a number of alignment results.
    Sorting criteria:
        1. aligned > un-aligned
        2. take higher alignment score
        3. take higher mapping quality
        4. if all the above are tied, pick one randomly

    Args:
        a list of alignment infos

    Returns:
        the best record from the provided list
    '''
    order = list(range(len(list_info)))
    list_is_unaligned = [info[0] != 1 for info in list_info]
    list_as = [info[0] for info in list_info]
    list_mapq = [info[1] for info in list_info]

    #: sort order: if_aligned > AS > MAPQ
    #: so perform sorting in reversed order
    list_mapq, list_as, list_is_unaligned, order = \
        zip(
            *sorted(zip(list_mapq, list_as, list_is_unaligned, order), reverse = True)
        )
    list_as, list_mapq, list_is_unaligned, order = \
        zip(
            *sorted(zip(list_as, list_mapq, list_is_unaligned, order), reverse = True)
        )
    list_is_unaligned, list_as, list_mapq, order = \
        zip(
            *sorted(zip(list_is_unaligned, list_as, list_mapq, order), reverse = True)
        )

    for i in range(1, len(list_info)):
        #: if not a tie
        if (list_is_unaligned[i] != list_is_unaligned[0]) or \
            (list_as[i] != list_as[0]) or \
            (list_mapq[i] != list_mapq[0]):
            order = order[:i]
            break
    return random.sample(order, 1)

def get_best_line(list_line):
    '''
    Given a number of SAM lines, select the best one.

    Args:
        a list of SAM lines
    
    Returns:
        the index of the best line, and the line itself
    '''
    list_info = []
    list_name = []
    for line in list_line:
        info = get_info_from_sam_line(line)
        list_name.append(info[0])
        try:
            assert info[0] == list_name[0]
        except:
            print (info)
            print (list_name)
            exit (1)
        list_info.append(info[1:])
    idx = compare_score_and_mapq(list_info)[0]
    return idx, list_line[idx]

def get_best_pair(list_line):
    '''
    Given a number of paired-end SAM records, select the best pair.

    Args:
        a list of SAM pairs. 
        `list_line[i]` and `list_line[len(list_line)/2 + i]` are in pair.
    
    Returns:
        the index of the best pair, and the pair itself
    '''
    list_info = []
    list_name = []
    for i, line in enumerate(list_line[: len(list_line)/2]):
        # info: QNAME, AS:i, MAPQ, flag
        info = get_info_from_sam_line(line, flag = True)
        info_mate = get_info_from_sam_line(list_line[i + len(list_line)/2], flag = True)

        # check if QNAMEs of a pair match and flags are reasonable
        try:
            assert info[0] == info_mate[0]
        except:
            print ('Error: read names between a pair do not match')
            print (info)
            print (info_mate)
        try:
            # flag 64 : first segment
            # flag 128 : second segment
            # A pair must have a first-seg read and a second-seg read
            assert ((info[3] & 64) ^ (info_mate[3] & 128)) or ((info[3] & 128) ^ (info_mate[3] & 64))
        except:
            print ('Error: segment information between a pair does not match')
            print (info)
            print (info_mate)

        list_name.append(info[0])

        # QNAMEs bewteen a set of matching SAM files must align
        try:
            assert info[0] == list_name[0]
        except:
            print ('Error: SAM records across input files are not aligned')
            print (info)
            print (list_name)
            exit (1)
        
        # compare sum of score and MAPQ
        pair_info = []
        pair_info.append(info[1] + info_mate[1])
        pair_info.append(info[2] + info_mate[2])
        list_info.append(pair_info)

    idx = compare_score_and_mapq(list_info)[0]
    return idx, list_line[idx], list_line[idx + len(list_line)/2]

def merge_core(list_f_sam, list_f_out, is_paired_end):
    '''
    This is the core function to perform merging.

    Args:
        list_f_sam: a list of SAM files to be merged
        list_f_out: 
            a list of output SAM files after merging. 
            `merge_core()` write results directly to the files
        is_paired_end: True if data is paired-end; False if single-end
    '''
    len_list = len(list_f_sam)
    list_read = []
    for line in list_f_sam[0]:
        # read all header lines for the first input SAM file
        if line[0] == '@':
            list_f_out[0].write(line)
            continue
        list_read.append(line)

        for i, f in enumerate(list_f_sam[1:]):
            f_line = f.readline()
            # read all header lines for the rest input SAM files
            while f_line[0] == '@':
                list_f_out[i+1].write(f_line)
                f_line = f.readline()
            list_read.append(f_line)
        
        try:
            # each read should be included by all SAM files
            if not is_paired_end:
                assert len(list_read) == len_list
            else:
                assert (len(list_read) == len_list) or (len(list_read) == 2 * len_list)
        except:
            print ('Error: number of alignments does not match')
            print ('len_read = {}, len_list = {}'.format(len(list_read), len_list))
            print (list_read)
            exit (1)

        if not is_paired_end:
            best_idx, best_line = get_best_line(list_read)
            list_f_out[best_idx].write(best_line)
            list_read = []
        else:
            if len(list_read) == 2 * len_list:
                best_idx, best_line1, best_line2 = get_best_pair(list_read)
                list_f_out[best_idx].write(best_line1)
                list_f_out[best_idx].write(best_line2)
                list_read = []
    return

def merge_incremental(args):
    '''
    Handles files I/O and processes arguments
    '''
    fn_sam = args.sam_list # 'to_merge.path'
    fn_ids = args.id_list # 'to_merge.id'
    fn_log = args.log # 'merged.path'
    prefix = args.prefix
    
    #: set random seed
    seed = args.rand_seed
    sys.stderr.write('Set random seed: {}\n'.format(seed))
    random.seed(seed)
    
    with open(fn_sam, 'r') as f:
        list_fn_sam = [] #: list of SAM names
        list_f_sam = [] #: list of opened SAM files
        for line in f:
            fn_sam = line.rstrip()
            list_fn_sam.append(fn_sam)
            list_f_sam.append(open(fn_sam, 'r'))
    with open(fn_ids, 'r') as f:
        list_ids = []
        for line in f:
            list_ids.append(line.rstrip())

    # number of files should match number of labels
    len_list = len(list_fn_sam)
    try:
        assert len_list == len(list_ids)
    except:
        print ('Error: numbers of files and labels do not match')
        exit (1)
    
    if fn_log != None:
        f_log = open(fn_log, 'w')
    list_f_out = []
    sys.stderr.write('output files: \n')
    for i in range(len(list_ids)):
        fn_out = prefix + '-' + list_ids[i] + '.sam'
        if fn_log != None:
            f_log.write(fn_out + '\n')
        sys.stderr.write(fn_out + '\n')
        list_f_out.append(open(fn_out, 'w'))

    merge_core(list_f_sam, list_f_out, args.paired_end)

    return

if __name__ == '__main__':
    args = parse_args()
    merge_incremental(args)
