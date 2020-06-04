# This is part of MiMiCPy

"""

This module contains the helper functions for the protein
module to strip water molecules within a certain radius

"""

import numpy as np
from ..parsers import pdb as hpdb

def _getSphere(coords):
    center = coords.mean(0)
    dist = np.linalg.norm(coords-center, axis=1)
    
    return center, np.amax(dist)

def _getCoords(pdb, chain=''):
    coords = []
    id_list = []
       
    for line in pdb.splitlines():
        vals = hpdb.readLine(line)
        if chain != '' and chain != vals['chainID']:
            continue
        
        idx = vals['serial']
        x = float(vals['x'])
        y = float(vals['y'])
        z = float(vals['z'])
        
        coords.append([x,y,z])
        id_list.append(idx)
        
    return np.array(coords), id_list

def coordsWithin(res, chain, hoh, dist):
    cen, rad = _getSphere(_getCoords(res, chain)[0])
    coords, idl = _getCoords(hoh)
    dist_vals = np.linalg.norm(coords - cen, axis=1)
    idx = np.where(dist_vals < (rad+dist))[0]
    return [idl[i] for i in idx]
	
