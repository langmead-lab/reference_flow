'''
This script reads population/superpopulation information from 
1000 Genomes Project and lists individuals wrt to population/
superpopulation.
'''

import argparse
import pandas as pd
# from plot_from_exps import read_ped
from utils import read_ped

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--fn_ped',
        help='ped file recoding population and family info'
    )
    parser.add_argument(
        '-sp', '--fn_superpopulation',
        help='superpopulation table'
    )
    parser.add_argument(
        '-op', '--out-prefix',
        help='prefix of output lists'
    )
    # parser.add_argument(
    #     '-c', '--cluster',
    #     help='txt file specifying unsupervised clusters'
    # )
    return parser.parse_args()

def process_db(
    fn_ped,
    fn_superpopulation
):
    '''
    Create dictionaries specifying mapping between indiv ID and
    1KG population/superpopulation labels.
    '''
    dict_indiv_pop = read_ped(fn_ped)
    
    df_superpop = pd.read_csv(
        fn_superpopulation,
        sep='\t',
        header=None,
        index_col=0
    )
    superpop_groups = df_superpop.groupby(2)
    
    dict_pop_superpop = {}
    for n, _ in superpop_groups:
        for i in superpop_groups.groups[n]:
            dict_pop_superpop[i] = n
    
    dict_indiv_superpop = {}
    for indiv in dict_indiv_pop.keys():
        dict_indiv_superpop[indiv] = dict_pop_superpop[dict_indiv_pop[indiv]]

    return dict_indiv_pop, dict_indiv_superpop

def list_indiv_from_pop(
    dict_indiv_pop,
    dict_indiv_superpop
):
    '''
    Reads dictionaries and returns dictionaries of lists
    where each list records indivs belong to a pop/superpop
    '''
    list_pop = list(set(dict_indiv_pop.values()))
    dict_list_pop = {}
    for pop in list_pop:
        list_pop = []
        for indiv in dict_indiv_pop.keys():
            if dict_indiv_pop[indiv] == pop:
                list_pop.append(indiv)
        dict_list_pop[pop] = list_pop
    
    list_superpop = list(set(dict_indiv_superpop.values()))
    dict_list_superpop = {}
    for superpop in list_superpop:
        list_superpop = []
        for indiv in dict_indiv_superpop.keys():
            if dict_indiv_superpop[indiv] == superpop:
                list_superpop.append(indiv)
        dict_list_superpop[superpop] = list_superpop
    return dict_list_pop, dict_list_superpop

def write_to_files(
    dict_list_pop,
    dict_list_superpop,
    out_prefix,
    num_indiv
):
    '''
    Writes dict of lists to files using given prefix
    '''
    count_pop_indiv = 0
    for pop in dict_list_pop.keys():
        fn_out = out_prefix + '_pop_' + pop + '.txt'
        f_out = open(fn_out, 'w')
        list_pop = dict_list_pop[pop]
        for indiv in list_pop:
            count_pop_indiv += 1
            f_out.write(indiv+'\n')
        f_out.close()
    assert count_pop_indiv == num_indiv
    
    count_superpop_indiv = 0
    for superpop in dict_list_superpop.keys():
        fn_out = out_prefix + '_superpop_' + superpop + '.txt'
        f_out = open(fn_out, 'w')
        list_superpop = dict_list_superpop[superpop]
        for indiv in list_superpop:
            count_superpop_indiv += 1
            f_out.write(indiv+'\n')
        f_out.close()
    assert count_superpop_indiv == num_indiv
    
    return

if __name__ == '__main__':
    args = parse_args()
    fn_ped = args.fn_ped
    fn_superpopulation = args.fn_superpopulation
    out_prefix = args.out_prefix

    dict_indiv_pop, dict_indiv_superpop = process_db(
        fn_ped,
        fn_superpopulation
    )
    assert len(dict_indiv_pop.keys()) == len(dict_indiv_superpop.keys())

    dict_list_pop, dict_list_superpop = list_indiv_from_pop(
        dict_indiv_pop,
        dict_indiv_superpop
    )

    write_to_files(
        dict_list_pop,
        dict_list_superpop,
        out_prefix,
        len(dict_indiv_pop.keys())
    )
