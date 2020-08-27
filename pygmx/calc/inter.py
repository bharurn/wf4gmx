import numpy as np
import pandas as pd
import xarray as xr
from tqdm.notebook import tqdm

def ids_dist(u, sele1, sele2, cutoff=3):
    mol = u.select_atoms(sele1)
    hyd = u.select_atoms(sele2)
    
    hpos = hyd.positions
    mpos = mol.positions
    pos = (mpos[:, np.newaxis] - hpos)
    
    reshaped = pos.reshape(-1, mpos.shape[1])
    dist = np.sum(np.abs(reshaped)**2,axis=1)**(1./2)
    dist = dist.reshape(pos.shape[0], pos.shape[1])
    
    ids = np.where(dist<cutoff)
    vals = dist[ids]
    
    conv = lambda atom: f"{atom.resid}{atom.resname}_{atom.name}"
    
    return list(map(conv, mol.atoms[ids[0]])), list(map(conv, hyd.atoms[ids[1]])), vals

def run(universe, sele1, sele2):
    _sele1 = []
    _sele2 = []
    time = []
    dist = np.empty(0)
    
    t0 = universe.trajectory[0].data['step']
    step = universe.trajectory[1].data['step'] - t0
    pbar = tqdm(total=len(universe.trajectory))
    
    universe.trajectory[0]
    for t in universe.trajectory:
        a,b,c = ids_dist(universe, sele1, sele2)
        _sele1 += a
        _sele2 += b
        time += [t.data['step']]*len(c)
        dist = np.append(dist, c)
        pbar.update(1)
    
    pbar.close()
    
    idx = pd.MultiIndex.from_arrays(arrays=[_sele1, _sele2, time], names=['Sele1','Sele2', "Time"])
    s = pd.Series(data=dist, index=idx)
    return xr.DataArray.from_series(s)
