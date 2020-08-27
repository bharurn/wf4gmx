# This is part of MiMiCPy

"""

This initializer script of mimicpy

"""

#from .system.ligand import NonStdLigand, StdLigand
#from .system.protein import Protein
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
    
gbl.host = Local('.', None)
gbl.logger = Logger(info=sys.stdout, debug=sys.stdout, debug2=None, warning=sys.stderr)
