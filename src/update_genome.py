'''
Updates a genome with a set of SNPs (and INDELS) and
write to a new file.
'''

import sys
import argparse
import random
import copy
from collections import OrderedDict

def get_mutation_type(orig, alts):
    '''
    Compare length of REF allele and ALT allele(s) and report variant type.

    Args:
        orig: ORIG genotype (string)
        alts: ALT genotype (string), if multi-allelic, split by ','
    Return:
        a string that can be the following:
            SNP: equal in length, length == 1
            MNP: equal in length, length > 1
            INDEL: not equal in length
            MULT: multiple ALT alleles
    Assertions:
        there is only one REF allele
    '''
    assert orig.count(',') == 0
    if alts.count(',') == 0:
        if len(orig) == len(alts) and len(orig) == 1:
            return 'SNP'
        elif len(orig) == len(alts) and len(orig) != 1:
            return 'MNP'
        elif len(orig) != len(alts):
            return 'INDEL'
    return 'MULT'
        
def get_allele_freq(info, num_haps, data_source, gnomad_ac_field):
    '''
    Returns allele frequency for a variant.
    Not using the "AF" attribute because it's calculated
    based on the entire 1KG population.
    '''
    attrs = info.split(';')
    #: for 1kg data, calculate allele frequency using phasing information
    if data_source == '1kg':
        for a in attrs:
            if a[:3] == 'AC=':
                try:
                    count = int(a[3:])
                #: when there are multiple alleles,
                #: use the highest frequency
                except:
                    a = a[3:].split(',')
                    inta = [int(i) for i in a]
                    count = max(inta)
        return float(count) / num_haps
    #: for genomad data, use pre-calculated allele frequency
    elif data_source == 'gnomad':
        for a in attrs:
            field = a.split('=')[0]
            if field == gnomad_ac_field:
                return float(a.split('=')[1]) / num_haps
    return -1


def process_vcf_header(line, indiv, f_vcf, data_source, is_ld):
    '''
    Process the header line of a VCF file

    Args:
        line (string): header line from the VCF file
        indiv (string): targeted sample (can be None)
        f_vcf (file): file that we are writing to
        data_source (string): project that generates the call set
        is_ld (boolean): if maintain local LD

    Returns:
        col (int/None): column index for `indiv`, None if `indiv` is not provided
        num_haps (int/None): number of samples in the call set, None if the call set not including phasing
        labels (list): split fields of `line`
    '''
    labels = line.rstrip().split('\t')
    # if `indiv` is set, select the corresponding column
    col = None
    num_haps = None
    if indiv != None:
        for i in range(9, len(labels)):
            if labels[i] == indiv:
                col = i
        if not col:
            print('Error! Couldn\'t find individual %s in VCF' % indiv)
            exit(1)
        f_vcf.write('\t'.join(labels[:9]) + f'\t{labels[col]}\n')
    else:
        if data_source == '1kg':
            if is_ld:
                f_vcf.write('\t'.join(labels[:9]) + '\tLD_SAMPLE\n')
            else:
                # skip sample genotype columns
                f_vcf.write('\t'.join(labels[:9]) + '\n')
        else:
            f_vcf.write(line)

    if data_source == '1kg':
        # calculate number of haplotypes (2 x number of samples)
        num_haps = 2 * (len(labels) - 9)

    return col, num_haps, labels


def write_to_fasta(dict_genome, out_prefix, suffix, line_width = 60):
    '''
    Write genome to a FASTA file
    '''
    f_fasta = open(f'{out_prefix}{suffix}.fa', 'w')

    # uncomment to show output in lexical order
    # for key in sorted(dict_genome.keys()):

    # show output following input order
    for key in dict_genome.keys():
        # write full contig name
        f_fasta.write(f'>{dict_genome[key][0]}\n')
        # write sequence, 60 chars per line
        for i in range(0, len(dict_genome[key][1]), line_width):
            f_fasta.write(''.join(dict_genome[key][1][i: i + line_width]) + '\n')
    f_fasta.close()


def update_allele(
    orig,
    alts,
    allele,
    indels,
    head,
    loc,
    f_var,
    hap,
    hap_str,
    offset,
    offset_other,
    chrom
):
    '''
    Update an allele
    '''

    if len(orig) != len(alts[allele-1]):
        v_type = 'INDEL'
    else:
        v_type = 'SNP'
    flag_skip = False
    if indels:
        # ignores conflicts or overlapped variants
        # but accepts overlapped INS
        if loc == head - 1 and (len(orig) < len(alts[allele - 1])):
            print ('Warning: overlapped INS at {0}, {1} for hap{2}'.format(loc, chrom, hap_str))
            new_offset = add_alt(hap, loc-1, orig, alts[allele-1], offset, True)
        elif loc >= head:
            new_offset = add_alt(hap, loc-1, orig, alts[allele-1], offset, False)
        else:
            flag_skip = True
            print ('Warning: conflict at {0}, {1} for hap{2}'.format(loc, chrom, hap_str))
    else:
        new_offset = 0
        hap[loc+offset-1] = alts[allele-1]

    if not flag_skip:
        f_var.write(
            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
            (hap_str, chrom, v_type, str(loc), str(loc+offset), orig, alts[allele-1], str(new_offset), str(offset_other) )
        )
        offset = new_offset
        head = loc + len(orig)
    
    return head, hap, offset


def update_variant(
    row,
    col,
    ld_hap,
    indiv,
    indels,
    data_source,
    is_ld,
    f_var,
    f_vcf,
    dict_genome,
    dict_genome_B,
    offsetA,
    offsetB,
    headA,
    headB,
    ld_indiv
):
    '''
    Update a variant, which may have one or two alleles
    '''

    chrom = row[0]
    loc = int(row[1])
    orig = row[3]
    alts = row[4].split(',')

    if is_ld:
        alleles = row[col].split('|')
        # always put the haplotype as alleleA for simplicity
        # the haplotype is not limited to gt[0] (it is randomly chose)
        alleleA = int(alleles[ld_hap])
        alleleB = 0
    elif indiv != None:
        # no `indiv` selected, take the allele
        alleles = row[col].split('|')
        alleleA = int(alleles[0])
        alleleB = int(alleles[1])
    else:
        #: always uses allele "1"
        alleleA = 1
        alleleB = 0
        
    if alleleA > 0:
        headA, dict_genome[chrom][1], offsetA = update_allele(
            orig=orig,
            alts=alts,
            allele=alleleA,
            indels=indels,
            head=headA,
            loc=loc,
            f_var=f_var,
            hap=dict_genome[chrom][1],
            hap_str='A',
            offset=offsetA,
            offset_other=offsetB,
            chrom=chrom
        )
    
    if alleleB > 0 and indiv != None:
        headB, dict_genome_B[chrom][1], offsetB = update_allele(
            orig=orig,
            alts=alts,
            allele=alleleB,
            indels=indels,
            head=headB,
            loc=loc,
            f_var=f_var,
            hap=dict_genome_B[chrom][1],
            hap_str='B',
            offset=offsetB,
            offset_other=offsetA,
            chrom=chrom
        )

    if (alleleA > 0) or \
        (alleleB > 0 and indiv != None):
            if data_source == '1kg':
                if indiv != None:
                    f_vcf.write('\t'.join(row[:9]) + f'\t{row[col]}\n')
                else:
                    if is_ld:
                        f_vcf.write('\t'.join(row[:9]) + f':SP\t{row[col]}:{ld_indiv}\n')
                    else:
                        f_vcf.write('\t'.join(row[:9]) + '\n')
            else:
                f_vcf.write(line)
    return headA, dict_genome, offsetA, headB, dict_genome_B, offsetB

def update_genome(
    indiv,
    dict_genome,
    vcf,
    out_prefix,
    indels,
    var_only,
    is_stochastic,
    block_size,
    is_ld,
    exclude_list,
    data_source,
    gnomad_ac_field,
    gnomad_pop_count,
    gnomad_af_th
):
    '''
    Handles variant updating for the entire genome. 
    This function mainly handles different updating settings and reads the VCF file.

    Variant updating follows the heirarchy:
        update_genome() -> update_variant() -> update_allele()
        - genome-level     - variant-level     - allele-level
    '''
    # assertions
    if is_ld:
        assert indiv == None
    # currently only supports 1000 Genomes ('1kg') and GnomAD ('gnomad') datasets
    assert data_source in ['1kg', 'gnomad']
    if data_source == 'gnomad':
        # GnomAD has no phasing information
        assert indiv == None
        assert is_ld == False
        assert exclude_list == ''

    '''
    ##fileformat=VCFv4.1
    '''
    if indiv != None:
        dict_genome_B = copy.deepcopy(dict_genome)
    else:
        dict_genome_B = None
    f_var = open(out_prefix + '.var', 'w')
    f_vcf = open(out_prefix + '.vcf', 'w')

    '''
    Format of a .var file:
    hap(A/B) chrm var_type ref_pos hap_pos offset
    '''
    if vcf != None:
        f = open(vcf, 'r')
    else:
        f = sys.stdin

    labels = None
    offsetA = 0
    offsetB = 0
    headA = 0
    headB = 0
    chrom = ''
    ld_hap = None
    ld_indiv = None

    if data_source == 'gnomad':
        num_haps = gnomad_pop_count
    else:
        num_haps = 0

    current_block_pos = 0
    rr = 0 # initial number for random.random()
    exclude_list = exclude_list.split(',')

    if is_ld:
        header_to_add_format = True
    else:
        header_to_add_format = False

    for line in f:
        #: Skip header lines
        if line[0] == '#' and line[1] == '#':
            if header_to_add_format and line.startswith('##FORMAT'):
                f_vcf.write('##FORMAT=<ID=SP,Number=1,Type=String,Description="Sample from which genotype is selected">\n')
            f_vcf.write(line)
            continue
        if line[0] == '#':
            col, num_haps, labels = process_vcf_header(line, indiv, f_vcf, data_source, is_ld)
            continue

        row = line.rstrip().split('\t')
        v_type = get_mutation_type(row[3], row[4])

        # switch to a new contig
        # we assume variants at different contigs are not interleaved
        if row[0] != chrom:
            headA = 0
            headB = 0
            offsetA = 0
            offsetB = 0
            chrom = row[0]
        loc = int(row[1])

        #: filter based on gnomad_af_th if it is set (gnomad only)
        if (not is_stochastic) and data_source == 'gnomad':
            freq = get_allele_freq(row[7], num_haps, data_source, gnomad_ac_field)
            if freq <= gnomad_af_th:
                continue

        #: no LD stochastic update for 1kg and gnomad
        if is_stochastic and is_ld == False:
            freq = get_allele_freq(row[7], num_haps, data_source, gnomad_ac_field)
            if freq < 0:
                continue
            #: only updates the random number when exceeding current block
            if loc >= current_block_pos + block_size:
                # print ('--update block--')
                # print ('prev rr = {0}, block_pos = {1}'.format(rr, current_block_pos))
                rr = random.random()
                current_block_pos = int(loc / block_size) * block_size
                # print ('updt rr = {0}, block_pos = {1}'.format(rr, current_block_pos))

            if rr > freq:
                # skip this allele
                continue
            # print ('selected, rr = {}'.format(rr), row[:2], freq)

        #: LD-preserving stochastic update for 1kg, this mode is not supported when using gnomad data
        if is_stochastic and is_ld and data_source == '1kg':
            if loc >= current_block_pos + block_size:
                while 1:
                    # randomly pick an individual
                    ld_indiv = random.choice(labels[9:])
                    # randomly pick a haplotype
                    ld_hap = random.choice([0,1])
                    if ld_indiv in exclude_list:
                        print ('exclude {0}: {1}-{2}'.format(current_block_pos, ld_indiv, ld_hap))
                        continue
                    current_block_pos = int(loc / block_size) * block_size
                    for i in range(9, len(labels)):
                        if labels[i] == ld_indiv:
                            col = i
                    if not col:
                        print('Error! Couldn\'t find individual %s in VCF' % indiv)
                        exit()    
                    break    

        if v_type == 'SNP' or (indels and v_type in ['INDEL', 'MULT']):
            headA, dict_genome, offsetA, headB, dict_genome_B, offsetB = update_variant(
                row = row,
                col = col,
                ld_hap = ld_hap,
                indiv = indiv,
                indels = indels,
                data_source = data_source,
                is_ld = is_ld,
                f_var = f_var,
                f_vcf = f_vcf,
                dict_genome = dict_genome,
                dict_genome_B = dict_genome_B,
                offsetA = offsetA,
                offsetB = offsetB,
                headA = headA,
                headB = headB,
                ld_indiv = ld_indiv
            )
            
    f_vcf.close()
    f_var.close()
    
    if not var_only:
        if indiv:
            # diploid output genome if `indiv` is set
            write_to_fasta(dict_genome, out_prefix, '_hapA')
            write_to_fasta(dict_genome_B, out_prefix, '_hapB')
        else:
            # haploid, no suffix
            write_to_fasta(dict_genome, out_prefix, '')

        # # write hapB sequence when `indiv` is set (diploid)
        # if indiv != None:
        #     fB = open(f'{out_prefix}_hapB.fa', 'w')
        #     # show output in lexical order
        #     # for key in sorted(dict_genome_B.keys()):
        #     for key in dict_genome_B.keys():
        #         # write full contig name
        #         fB.write(f'>{dict_genome_B[key][0]}\n')
        #         # write sequence, 60 chars per line
        #         for i in range(0, len(dict_genome_B[key][1]), 60):
        #             fB.write(''.join(dict_genome_B[key][1][i: i+60]) + '\n')
        #     fB.close()
    f.close()

def add_alt(genome, loc, orig, alt, offset, overlap_ins):
    '''
    loc here is the index for str
    i.e., loc = vcf_loc -1
    '''
    loc += offset

    if len(orig) == 1 and len(alt) == 1:
        # SNP
        # if genome[loc] != orig:
        #    print (loc, genome[loc], alt)
        genome[loc] = alt
    elif len(orig) > len(alt):
        # Deletion
        for i in range(len(alt)):
            genome[loc+i] = alt[i]
        del genome[loc+len(alt):loc+len(orig)]
        offset -= (len(orig) - len(alt))
    elif len(alt) > len(orig):
        #: Insertion
        
        # if overlap_ins:
        #     for i in range(len(orig)):
        #         if genome[loc+i] != alt[i]:
        #             print ('Warning: genome ({0}) differs from ALT ({1}) at {2}'.format(genome[loc+i], alt[i], loc))
        
        #: don't replace if overlap_ins is True
        if not overlap_ins:
            for i in range(len(orig)):
                genome[loc+i] = alt[i]
        genome[loc+len(orig):loc+len(orig)] = list(alt[len(orig):])
        offset += len(alt) - len(orig)
    else:
        # len(orig)=len(alt)>1 : Weird combination of SNP/In/Del
        for i in range(len(alt)):
            genome[loc+i] = alt[i]

    return offset

def read_genome(fn_genome):
    '''
    Read a fasta file and returns a dictionary storing contig sequences

    Args:
        fn_genome: file name of the genome (in fasta format)

    Return:
        dict_genome (dict)
            key: short contig name (from '>' to the first space)
            value: [full contig name (string), sequence (list)]
                sequence is represented in lists for easier modification later
    '''
    if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
        # Use OrderedDict to maintain dictionary oder as input
        dict_genome = OrderedDict()
    else:
        # Python 3.6's dict is now ordered by insertion order
        dict_genome = {}
    
    f_genome = open(fn_genome, 'r')
    curr_contig = ''
    curr_seq = ''
    full_contig_name = ''
    for line in f_genome:
        line = line.rstrip()
        if line[0] == '>':
            # update previous contig
            if curr_contig != '':
                dict_genome[curr_contig] = [full_contig_name, list(curr_seq[:])]
            full_contig_name = line[1:]
            curr_contig = line.split()[0][1:]
            curr_seq = ''
        else:
            curr_seq += line
    if curr_contig != '':
        dict_genome[curr_contig] = [full_contig_name, list(curr_seq[:])]

    return dict_genome

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-r', '--ref', type=str, required=True, help='Path to fasta file containing reference genome'
    )
    parser.add_argument(
        '-v', '--vcf', type=str, help="Path to VCF file containing mutation information"
    )
    # I'd like to deprecate this
    parser.add_argument(
        '-c', '--chrom', type=str, help="Chromosome to process (deprecated)"
        # '-c', '--chrom', type=str, required=True, help="Chromosome to process"
    )
    parser.add_argument(
        '-op', '--out-prefix', type=str, required=True, help="Path to output prefix"
    )
    parser.add_argument(
        '-s', '--name', type=str, help="Name of individual in VCF to process; leave blank to allow all variants [None]"
    )
    parser.add_argument(
        '-i', '--include-indels', action='store_true', help="Set to extract both SNPs and INDELs [Off]"
    )
    parser.add_argument(
        '-S', '--stochastic', action='store_true', help="Set to enable stochastic flipping [Off]"
    )
    parser.add_argument(
        '-rs', '--rand-seed', help="random seed for controlled randomness [None]"
    )
    parser.add_argument(
        '--var-only', action='store_true', help="Set to report .var file only (no .fa output) [Off]"
    )
    parser.add_argument(
        '-b', '--block-size', type=int, default=1, help="Size of block for stochastic update [1]"
    )
    parser.add_argument(
        '-l', '--ld', action='store_true', help="Set to enable pseudo-LD blocking [Off]"
    )
    parser.add_argument(
        '-ex', '--exclude-name', type=str, default='', help="Name of individuals in VCF to exclude; separate by comma ['']"
    )
    parser.add_argument(
        '-d', '--data-source', type=str, default='1kg', help="Source of population genomic data, currently support '1kg' and 'gnomad' ['1kg']"
    )
    parser.add_argument(
        '--gnomad-ac-field', type=str, default='AF',
        help="GnomAD allele count field; activated only in stochastic mode; \
            can be changed depending on popultion of interest ['AC']"
    )
    parser.add_argument(
        '--gnomad-pop-count', type=int, help="Size of GnomAD population [INT]"
    )
    parser.add_argument(
        '--gnomad-af-th', type=float, default=0,
        help="GnomAD allele frequency threshold. Variants with frequency lower than this value \
            will not be updated [0]"
    )

    args = parser.parse_args(sys.argv[1:])

    if args.name == None:
        print ('Note: no individual specified, all variants are considered')
    # if args.stochastic == 1:
    if args.stochastic:
        print ('Note: stochastic update is enabled')
        if args.rand_seed:
            random.seed(args.rand_seed)
            print ('Set random seed: {}'.format(args.rand_seed))

    dict_genome = read_genome(args.ref)
    update_genome(
        indiv = args.name,
        dict_genome = dict_genome,
        vcf = args.vcf,
        out_prefix = args.out_prefix,
        indels = args.include_indels,
        var_only = args.var_only,
        is_stochastic = args.stochastic,
        block_size = args.block_size,
        is_ld = args.ld,
        exclude_list = args.exclude_name,
        data_source = args.data_source,
        gnomad_ac_field = args.gnomad_ac_field,
        gnomad_pop_count = args.gnomad_pop_count,
        gnomad_af_th = args.gnomad_af_th
    )
