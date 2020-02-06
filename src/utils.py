import pandas as pd

def read_ped(ped_fn):
    '''
    Reads the ped file and returns a dict where keys are indivs and 
    values are corresponding populations
    '''
    ped_df = pd.read_csv(ped_fn, sep='\t')
    dd = ped_df[['Individual ID', 'Population']]

    popd = {}
    for i in range(dd.shape[0]):
        key = dd['Individual ID'][i]
        value = dd['Population'][i]
        popd[key] = value

    return popd
