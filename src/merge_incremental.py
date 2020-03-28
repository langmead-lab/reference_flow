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
    args = parser.parse_args()
    return args

def get_name_as_mapq(line):
    '''
    Read one line from SAM and return QNAME, alignment score (AS:i) and MAPQ

    Args:
        a line from SAM

    Returns:
        QNAME, AS:i, MAPQ
    '''
    if line[0] == '@':
        return 'header', None
    line = line.split()
    name = line[0]
    mapq = int(line[4])
    score = 1 #: represents unmapped
    for i in line:
        if i.startswith('AS:i'):
            score = int(i.split(':')[-1])
    return name, score, mapq

def read_line(istream):
    while True:
        line = istream.readline()
        if line and line[0] != '@':
            return line

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
        info = get_name_as_mapq(line)
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

if __name__ == '__main__':
    args = parse_args()
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
    
    for line in list_f_sam[0]:
        if line[0] == '@':
            list_f_out[0].write(line)
            continue
        list_read = [line]
        for i, f in enumerate(list_f_sam[1:]):
            f_line = f.readline()
            while f_line[0] == '@':
                list_f_out[i+1].write(f_line)
                f_line = f.readline()
            # f_line = read_line(f, i)
            list_read.append(f_line)

        # each read should be included by all SAM files
        try:
            assert len(list_read) == len_list
        except:
            print ('Error: number of alignments does not match')
            print (list_read)
            exit (1)
        best_idx, best_line = get_best_line(list_read)
        list_f_out[best_idx].write(best_line)
    
