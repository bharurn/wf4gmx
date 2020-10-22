# This is part of MiMiCPy

"""

This initializer script of mimicpy

"""

from .system.ligand import NonStdLigand, StdLigand
from .system.protein import Protein
from .handles.prepare import Prepare
from .handles.analysis import Analyze
from .handles.simulate import MD, MiMiC
from .scripts.slurm import Slurm
from .scripts.mdp import MDP
from .shell.shell import Local, Remote
from ._global import _Global as gbl
from .utils.logger import Logger
from .utils.errors import MiMiCPyError
from . import calc
import sys

def setHost(dirc='.', *args, path=None):
    """Wrapper function to set-up host"""

    closeHost()
    
    if ':' not in dirc:
        gbl.host = Local(dirc, path, *args)
    else:
        gbl.host = Remote(dirc, path, *args)

def getHost():
    """Convenience function to return host"""
    return gbl.host

def setEnv(**kwargs):
    """Set all executable paths/names in _Global"""
    for k,v in kwargs.items():
        if hasattr(gbl, k) and (k != 'host' or k != 'logger'):
            setattr(gbl, k, v)
        else:
            raise MiMiCPyError(f"{k} is not an enviornment executable/path!")

def setLogger(level, redirect=sys.stdout):
    """Set the logger level and stream"""
    if level == 0: # no output
        gbl.logger.info = None
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 1: # only info
        gbl.logger.info = redirect
        gbl.logger.debug = None
        gbl.logger.debug2 = None
    elif level == 2: # some debug commands
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = None
    elif level == 3: # verbose debug commands
        gbl.logger.info = redirect
        gbl.logger.debug = redirect
        gbl.logger.debug2 = redirect
        
def redirectWarnings(redirect):
    gbl.logger.warning = redirect
        
def closeHost():
    """Convenience function to close host"""
    gbl.host.__del__()
    
def fromRCSB(pdbid, chains=None, asRaw=False):
   import urllib.request as req
   from collections import defaultdict
   from .. import _hndlpdb as hpdb
   
   gbl.logger.write('info', f"Accessing PDB {pdbid} from RCSB database..")
   gbl.logger.write('debug', "Downloading PDB file..")
      
   url = f'http://files.rcsb.org/view/{pdbid}.pdb'
   gbl.logger.write('debug', f"Openeing {url}..")
   response = req.urlopen(url)
   #response.raise_for_status()
   
   pdb = response.read().decode('utf-8')
   
   ligs = defaultdict(str)
   
   gbl.logger.write('info', "Downloading ligands..")
   
   for l in pdb.splitlines():
        vals = hpdb.readLine(l)
        
        if vals['record'].strip() == 'HET':
            
            if chains is not None and vals['content'].split()[1] not in chains:
                continue
            
            v = vals['content'].split()[:3]
    
            query = '_'.join(v)
            
            gbl.logger.write('debug', f"Downloading ligands {v[0]}, chain {v[1]}")
    
            url = f"https://files.rcsb.org/cci/download/{pdbid}_{query}_NO_H.sdf"
            gbl.logger.write('debug', f"Openeing {url}..")
            resp = req.urlopen(url)
            #resp.raise_for_status()

            ligs[v[0]] += resp.read().decode('utf-8')
            
        elif vals['record'] == 'HETNAM': break
        
   #pdb_ = ''
   
   #for line in pdb.splitlines():
   #     vals = hpdb.readLine(line)
        
    #    if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
    #        if chains is None or vals['chainID'] in chains:
    #            pdb_ += line + '\n'
   
   if asRaw:
       gbl.logger.write('info', "Returned as raw data..")
       return pdb, ligs
   else:
       with open(f"{pdbid}.pdb", 'w') as f:
           f.write(pdb)
       for k,v in ligs.items():
           with open(f"{k}.sdf", 'w') as f: f.write(v)
           
       gbl.logger.write('info', "Saved to files..")
    
gbl.host = Local('.', None)
gbl.logger = Logger(info=sys.stdout, debug=sys.stdout, debug2=None, warning=sys.stderr)
