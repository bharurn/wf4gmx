# This is part of MiMiCPy

"""

This module contains the _Global which holds all global variables

"""

class _Global:
    """
    Class to hold all module scope varibale
    All variables are static members, so changes in their values are reflected everywhere
    """
    host = None
    logger = None
    gmx = None
    cpmd = 'cpmd.x'
    cpmd_pp = None
    gauss = 'g09'
    obabel = 'obabel'
    antechamber = 'antechamber'
    parmchk = 'parmchk2'
    tleap = 'tleap'
    acpype = 'acpype'