# This is part of MiMiCPy

"""

This module contains the helper functions for the ligand
module to generate topology and position restratint information
using the Generalized Amber Force Field (GAFF)

"""

from mimicpy._global import _Global as _global
from mimicpy.utils.errors import ExecutionError

def _cleanprep(mol, prep_to_pdb):
    prep = _global.host.read(f'{mol}.prep')
    
    for i in prep_to_pdb:
        prep = prep.replace(i, prep_to_pdb[i])

    return prep

def _cleanItp(itp, mol):
    atoms = False
    atomtypes = False
    notwrite = False
    itp_str = ""

    mol_ = f"{mol}_" # MOL_ string
    
    hvy = []
    
    atm_types = [] # store atom types for MPR
    
    for line in itp.splitlines()[1:]:
    
        if "[ atomtypes ]" in line:
            atomtypes = True # start atomtypes
        elif atomtypes:
            if line.strip() == '': # end atomtypes
                atomtypes = False
            elif not line.strip()[0] == ';': # if line is not comment
                splt = line.split() # split line
                #and stich it back together with splt[0] and splt[1], which are the atom types
                #repalced by mol_atom-type
                atm_types.append(mol_ + splt[0])
                line = ' '+ mol_ + splt[0] + ' '*7 + mol_ + splt[1] + line[12:]
            
        if "[ atoms ]" in line:
            atoms = True # start atoms
        elif atoms:
            if line.strip() == '' or line == '[ bonds ]':
                atoms = False # end atoms
            elif not line.strip()[0] == ';' :
                splt = line.split()
                #same as above
                line = line[:9] + mol_ + splt[1] + line[11:]
                
                if 'H' not in splt[4].upper():
                    hvy.append(splt[0])
        
        if "[ defaults ]" in line or "[ system ]" in line or "[ molecules ]" in line:
            notwrite = True
        
        if notwrite and line.strip() == '': notwrite = False
        
        if notwrite == False:
            itp_str += line+'\n' # collate line into itp_string
 
    return itp_str, hvy, atm_types

def _getposre(hvy):
    posre = "[ position_restraints ]"
    posre = "[ position_restraints ]\n" + '\n'.join([f"  {i}   1  1000 1000 1000" for i in hvy])
    return posre

def do(mol, conv):
    _global.logger.write('info', "Generating topology for the ligand..")
    
    prep = f"{mol}.prep"
    
    _global.host.checkFile(prep)
    
    if conv != {}:
        _global.logger.write('info', f"Mapping prep of {mol}.prep atom names to that of pdb using dictionary provided..")
        
        _global.host.write(_cleanprep(mol, conv), "params.prep")
        
        _global.logger.write('debug', "Ouput saved to params.prep")
        
        prep = "params.prep"
        
    _global.logger.write('info', f"Running AmberTools parmchk on {prep}..")
    
    log = _global.host.run(f'{_global.parmchk} -i {prep} -f prepi -o params.frcmod..')
    
    _global.logger.write('debug2', log)
    
    _global.logger.write('debug', "Output saved to params.frcmod..")
    
    tleap_in = "source leaprc.gaff\n"
    tleap_in += f"loadamberprep {prep}\n"
    #tleap_in += f"loadamberparams params.frcmod\n"
    
    if _global.host.checkFile(f"{mol}.frcmod", throw=False):
        tleap_in += f"loadamberparams {mol}.frcmod\n"
    
    tleap_in += f"saveamberparm {mol} {mol}.prmtop {mol}.inpcrd\nquit"
    
    _global.logger.write('info', "Running AmberTools LEaP..")
    cmd = f'tleap -f -'
    output = _global.host.run(cmd, stdin=tleap_in)
    
    log = '\n Running'+ tleap_in + '..\n' + output + '\n'
    _global.logger.write('debug2', log)
    
    if _global.host.fileExists(f'{mol}.prmtop') == False or _global.host.fileExists(f'{mol}.inpcrd') == False:
        raise ExecutionError(cmd, output)
    
    _global.logger.write('debug', f"Output saved to {mol}.prmtop and {mol}.inpcrd..")
    
    _global.logger.write('info', "Converting to Gromacs topology using Acpype..")
    cmd = f'acpype -p {mol}.prmtop -x {mol}.inpcrd'
    output = _global.host.run(cmd)
    
    log = f"Running {cmd}..\n" + output + '\n'
    _global.logger.write('debug2', log)
    
    if _global.host.fileExists(f'{mol}_GMX.gro') == False or _global.host.fileExists(f'{mol}_GMX.top') == False:
        raise ExecutionError(cmd, output)
        
    _global.logger.write('debug', f"Output saved to {mol}_GMX.gro and {mol}_GMX.top..")
    
    _global.logger.write('info', "Cleaning Acpype output..")
    itp = _global.host.read(f'{mol}_GMX.top')
    
    _global.logger.write('debug', "Renaming atoms to avoid conflict..")
    itp, hvy, atm_types = _cleanItp(itp, mol) # rename atoms, remove uneeded sections, get list of heavy atoms for posre
    
    _global.logger.write('debug', "Writing position restraint data..")
    posre = _getposre(hvy)
    
    itp += f'\n; Include Position restraint file\n#ifdef POSRES\n#include "posre_{mol}.itp"\n#endif'
    
    return itp, posre, atm_types

    
