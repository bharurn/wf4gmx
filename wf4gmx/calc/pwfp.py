import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
from .inter import ids_dist
    
def pairwise_fingerprint(u, lig, receptor='protein', nframes=None, nres=3, cutoff=10):    
    u.trajectory[0]
    lst = []
    if nframes is None: nframes = len(u.trajectory)
    pbar = tqdm(total=nframes)
    for i in u.trajectory[:nframes]:
        da1 = ids_dist(u, lig, receptor, cutoff=cutoff)
        df = pd.DataFrame(da1).T
        df.columns = ['Ligand', 'Receptor', 'Dist']
        
        ## find min dist between each lig atom and each receptor
        # get min dist for each lig receptro atom combo
        df = df.groupby(['Ligand', 'Receptor']).min()
        # drop receptor column
        df = df.reset_index().drop('Receptor', axis=1)
        
        ## find nres min values for each lig atom
        # func to sort series and return first nres min values
        # rest of the series values are but as nan, to maintain size
        min_n = lambda x: x.sort_values()[:nres].append(pd.Series([np.nan]*(len(x)-nres)))
        
        # call func for each ligand atom, and drop nan vals
        vals = df.groupby('Ligand').transform(min_n).dropna()['Dist'].to_numpy()
        
        lst.append(vals)
        pbar.update(1)
    
    pbar.close()

    return np.stack(lst)

def cluster(fprint, max_d=2.5, trunc=60):
    from scipy.cluster import hierarchy
    import matplotlib.pyplot as plt

    linked = hierarchy.linkage(fprint, 'single')

    labelList = range(0, fprint.shape[0])

    plt.figure(figsize=(10, 7))
    hierarchy.dendrogram(linked,
                orientation='top',
                truncate_mode='lastp',
                p=trunc,
                labels=labelList,
                distance_sort='descending',
                show_leaf_counts=True)
    plt.axhline(y=max_d, c='k')
    plt.show()
    
    clus = hierarchy.fcluster(linked, max_d, criterion='distance')
    cluster = pd.DataFrame([clus, labelList]).T
    cluster.columns=['cluster', 'frame']
    return cluster