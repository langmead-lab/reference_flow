'''
Updates a genome with a set of SNPs (and INDELS) and
write to a new file.
'''

import sys
import argparse
import random

def get_mutation_type(orig, alts):
    '''
    Compare length of REF allele and ALT allele(s) and report variant type:
        SNP: equal in length, length == 1
        MNP: equal in length, length > 1
        INDEL: not equal in length
        MULT: multiple ALT alleles
        (assertion fails if more than one REF allele)
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
        
# def get_mutation_type(info):
#     '''
#     Returns the value of the VT attribute from the info field
#     '''

#     attrs = info.split(';')
#     for a in attrs:
#         if a[:3] == 'VT=':
#             return a[3:]

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
    if len(orig) != len(alts[allele-1]):
        type = 'INDEL'
    else:
        type = 'SNP'
    flag_skip = False
    if indels:
        #: ignores conflicts or overlapped variants
        #: but accepts overlapped INS
        if loc == head-1 and (len(orig) < len(alts[allele-1])):
            print ('Warning: overlapped INS at {0} for hap{1}'.format(loc, hap_str))
            new_offset = add_alt(hap, loc-1, orig, alts[allele-1], offset, True)
        elif loc >= head:
            new_offset = add_alt(hap, loc-1, orig, alts[allele-1], offset, False)
        else:
            flag_skip = True
            print ('Warning: conflict at {0} for hap{1}'.format(loc, hap_str))
    else:
        new_offset = 0
        hap[loc+offset-1] = alts[allele-1]

    if not flag_skip:
        f_var.write(
            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
            (hap_str, chrom, type, str(loc), str(loc+offset), orig, alts[allele-1], str(new_offset), str(offset_other) )
        )
        offset = new_offset
        head = loc + len(orig)
    
    return head, hap, offset

def update_genome(
    indiv,
    seq,
    label,
    vcf,
    chrom,
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
    #: assertions
    if is_ld:
        assert indiv == None

    #: currently only supports 1000 Genomes ('1kg') and GnomAD ('gnomad') datasets
    assert data_source in ['1kg', 'gnomad']
    if data_source == 'gnomad':
        #: GnomAD has no phasing information
        assert indiv == None
        assert is_ld == False
        assert exclude_list == ''

    '''
    ##fileformat=VCFv4.1
    '''
    hapA = list(seq[:])
    hapB = list(seq[:])
    # if var_only == 0:
    if not var_only:
        if indiv != None:
            fA = open(out_prefix + '_hapA.fa', 'w')
            split_label = label.split()
            label_A = split_label[0] + 'A ' + ' '.join(split_label[1:]) + '\n'
            fA.write(label_A)
            fB = open(out_prefix + '_hapB.fa', 'w')
            label_B = split_label[0] + 'B ' + ' '.join(split_label[1:]) + '\n'
            fB.write(label_B)
        else:
            fA = open(out_prefix + '.fa', 'w')
            fA.write(label)

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
    line_id = 0
    offsetA = 0
    offsetB = 0
    headA = 0
    headB = 0
    if data_source == 'gnomad':
        num_haps = gnomad_pop_count
    else:
        num_haps = 0

    current_block_pos = 0
    rr = 0 # initial number for random.random()
    exclude_list = exclude_list.split(',')

    for line in f:
        #: Skip header lines
        if line[0] == '#' and line[1] == '#':
            f_vcf.write(line)
            continue
        if line[0] == '#':
            f_vcf.write(line)
            labels = line.rstrip().split('\t')
            if data_source == '1kg':
                num_haps = 2 * (len(labels) - 9)
            #: if "indiv" is set, select corresponding columns
            if indiv != None:
                col = None
                for i in range(9, len(labels)):
                    if labels[i] == indiv:
                        col = i
                if not col:
                    print('Error! Couldn\'t find individual %s in VCF' % indiv)
                    exit(1)
            continue
        row = line.rstrip().split('\t')
        # type = get_mutation_type(row[7])
        type = get_mutation_type(row[3], row[4])

        if row[0] != chrom and row[0] != 'chr' + chrom:
            continue
        loc = int(row[1])

        #: filter based on gnomad_af_th if it is set (gnomad only)
        if (not is_stochastic) and data_source == 'gnomad':
            freq = get_allele_freq(row[7], num_haps, data_source, gnomad_ac_field)
            if freq <= gnomad_af_th:
                continue

        #: no LD stochastic update for 1kg and gnomad
        if is_stochastic and is_ld == False:
            freq = get_allele_freq(row[7], num_haps, data_source, gnomad_ac_field)
            # print ('freq', freq)
            if freq < 0:
                # print ('Warning! gnomad_ac_field ({}) is not found'.format(gnomad_ac_field))
                # print (line)
                continue
            #: only updates the random number when exceeding current block
            if loc >= current_block_pos + block_size:
                # print ('--update block--')
                # print ('prev rr = {0}, block_pos = {1}'.format(rr, current_block_pos))
                rr = random.random()
                current_block_pos = int(loc / block_size) * block_size
                # print ('updt rr = {0}, block_pos = {1}'.format(rr, current_block_pos))

            if rr > freq:
                continue
            # print ('selected, rr = {}'.format(rr), row[:2], freq)

        #: LD-preserving stochastic update for 1kg, this mode is not supported when using gnomad data
        if is_stochastic and is_ld and data_source == '1kg':
            if loc >= current_block_pos + block_size:
                while 1:
                    ld_indiv = random.choice(labels[9:])
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

        #: supports tri-allelic
        # if type == 'SNP' or (indels and type in ['INDEL', 'SNP,INDEL']):
        if type == 'SNP' or (indels and type in ['INDEL', 'MULT']):
            orig = row[3]
            alts = row[4].split(',')

            if is_ld:
                alleles = row[col].split('|')
                alleleA = int(alleles[ld_hap])
                alleleB = 0
            elif indiv != None:
                alleles = row[col].split('|')
                alleleA = int(alleles[0])
                alleleB = int(alleles[1])
            else:
                #: always uses allele "1"
                alleleA = 1
                alleleB = 0
                
            if alleleA > 0:
                headA, hapA, offsetA = update_allele(
                    orig=orig,
                    alts=alts,
                    allele=alleleA,
                    indels=indels,
                    head=headA,
                    loc=loc,
                    f_var=f_var,
                    hap=hapA,
                    hap_str='A',
                    offset=offsetA,
                    offset_other=offsetB,
                    chrom=chrom
                )
            
            if alleleB > 0 and indiv != None:
                headB, hapB, offsetB = update_allele(
                    orig=orig,
                    alts=alts,
                    allele=alleleB,
                    indels=indels,
                    head=headB,
                    loc=loc,
                    f_var=f_var,
                    hap=hapB,
                    hap_str='B',
                    offset=offsetB,
                    offset_other=offsetA,
                    chrom=chrom
                )

            if (alleleA > 0) or \
                (alleleB > 0 and indiv != None):
                f_vcf.write(line)
            line_id += 1
    
    f_vcf.close()
    f_var.close()
    
    if not var_only:
        for i in range(0, len(hapA), 60):
            fA.write(''.join(hapA[i:i+60])  + '\n')
        fA.close()

        if indiv != None:
            for i in range(0, len(hapB), 60):
                fB.write(''.join(hapB[i:i+60])  + '\n') 
            fB.close()
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

def read_chrom(ref, chrom):
    with open(ref, 'r') as f_in:
        label = None
        seq = ''
        for line in f_in:
            if line[0] == '>':
                if label:
                    return label, seq
                curr_chrom = line[1:].split(' ')[0]
                curr_chrom = curr_chrom.rstrip()
                if (curr_chrom == chrom) or (curr_chrom == 'chr' + chrom):
                    label = line
            elif label:
                seq += line.rstrip()
        if label == None:
            print ('Error: no matched chromosome. Label = {}'.format(chrom))
            exit ()
        return label, seq

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
    parser.add_argument(
        '-c', '--chrom', type=str, required=True, help="Chromosome to process"
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
        print ('Note: no individual specified, all variants in chrom %s are included' % args.chrom)
    # if args.stochastic == 1:
    if args.stochastic:
        print ('Note: stochastic update is enabled')
        if args.rand_seed:
            random.seed(args.rand_seed)
            print ('Set random seed: {}'.format(args.rand_seed))

    label, genome = read_chrom(args.ref, args.chrom)
    update_genome(
        indiv = args.name,
        seq = genome,
        label = label,
        vcf = args.vcf,
        chrom = args.chrom,
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
