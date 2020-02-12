def build_dict_indiv_to_pop(fn):
    dict_indiv_to_pop = {}
    df = pd.read_csv(fn, sep='\t')
    list_sample = df['Individual ID']
    list_pop = df['Population']
    assert len(list_sample) == len(list_pop)
    assert len(list_sample) == len(set(list_sample))
    for i, s in enumerate(list_sample):
        dict_indiv_to_pop[s] = list_pop[i]
    return dict_indiv_to_pop

def build_dict_pop_to_spop(fn):
    dict_pop_to_spop = {}
    df = pd.read_csv(fn, sep='\t', header=None)
    list_spop = df[2]
    list_pop = df[0]
    for i, p in enumerate(list_pop):
        dict_pop_to_spop[p] = list_spop[i]
    return dict_pop_to_spop
