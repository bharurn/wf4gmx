# This is part of MiMiCPy

"""

This module contains the helper functions for the ligand
module to protonate ligands and convert sdf to pdf

"""

from . import _hndlpdb as hpdb
from .._global import _Global as _global
import re
from mimicpy.utils.errors import ExecutionError

def _cleanPdb(sdf, pdb):
    _global.logger.write('debug', "Assigning correct atom names..")
    
    #####Read sdf file#####
    l_cnt = 0
    start = False
    no = -1
    resname = ""
    atm_names = [] # store atoms names
    elems = [] # store element list for MPT
    for line in sdf.splitlines():
        if l_cnt == 3:
            no = int(line.split()[0])
        elif no != -1 and l_cnt < no+3: # read coordinates
            elems.append(line.split()[3]) #4th value is the element symbol
        elif "A    1" in line: # start reading atom names
            start = True
            val = line.split()
            if len(val) == 1:
                # get atom names to be added to the pdb file
                atm_names.append(val[0])
        elif "M  END" in line:
            start = False # stop reading atom names
        elif start:
            val = line.split()
            if len(val) == 1:
                # get atom names to be added to the pdb file
                atm_names.append(val[0])
        # search for ChemCompId tag to get ligand name
        elif "ChemCompId" in line:
            resname = "start"
        elif resname == "start":
            resname = line.strip()
        elif '$$$$' in line:
            l_cnt = -1 # reset value to zero for next chain
            no = -1
        l_cnt += 1
            

    #####END#####
    
    pdb_str = ""
    atm_i = 0 # counter for atom name
    h_i = 1 # counter for new hydrogen atom names
    # hydrogen atom names are added as H1, H2, H3, ....
    for line in pdb.splitlines():
        header = line.split()[0]
        if header == "HETATM":
            line = hpdb.editLine(line, name=atm_names[atm_i], resName=resname)
            atm_i += 1
        elif header == "ATOM":
        #ATOM is hydrogen atom, so add H1, H2...
            hstr = "H" + str(h_i) # create hydrogen atom name
            line = hpdb.editLine(line, name=hstr, resName=resname) # add the res name also
            h_i += 1
        
        pdb_str += line + '\n'

    return pdb_str, resname, elems



def _multChains(pdb):
    pdb_str = ""
    n_chains = 0
    stack = []
    stack_idx = 0
    
    _global.logger.write('debug', "Assigning correct chain IDs..")
    
    for line in pdb.splitlines():
        splt = hpdb.readLine(line)
        
        if splt['record'] == "COMPND":
            n_chains += 1
            stack_idx = 0
            
        if splt['record'] == "HETATM" or splt['record'] == "ATOM":
            
            line = hpdb.editLine(line, chainID=chr(ord('A')+n_chains-1), resSeq=str(n_chains))
            if splt['element'] == 'H':
                if n_chains == 1:
                    stack.append(hpdb.readLine(line)['name'])
                else:
                    line = hpdb.editLine(line, name=stack[stack_idx])
                    stack_idx += 1
            
        pdb_str += line + '\n'

    return pdb_str, n_chains

def _common(sdf, pdb):
    pdb, resname, elems = _cleanPdb(sdf, pdb)
    pdb, chains = _multChains(pdb)
    pdb = '\n'.join(re.findall(r"^(ATOM.*|HETATM.*)", pdb, re.MULTILINE))
    return pdb, chains, resname, elems

def _ob(cmd, stdin):
    pdb = _global.host.run(cmd, stdin=stdin)
    reg = re.compile(r"^(COMPND (.*\n)*)", re.MULTILINE) # separate out and err
    
    if reg == None:
        raise ExecutionError(cmd, reg)
    
    out = reg.search(pdb).groups()[0]
    err = pdb.replace(out, '')
    
    _global.logger.write('debug2', f"Running {cmd}..\n"+err)
    
    return out

def clean_sdf(sdf):
    _global.logger.write('debug', f"Cleaning sdf..")
    return re.sub(r'/A    1/(\n.*?)*?/M  END/', '', sdf, flags=re.MULTILINE)
    
def do(sdf, pH):
    sdf_text = clean_sdf(sdf)
    _global.logger.write('info', "Converting to pdb using Openbabel")
    pdb = _ob(f'{_global.obabel} -isdf -opdb', sdf_text)
    _global.logger.write('info', 'Adding hydrogens at pH={pH}')
    pdb = _ob(f'{_global.obabel} -ipdb -p {pH} -opdb', stdin=pdb)
    
    return _common(sdf, pdb)

def donoH(sdf):
    sdf_text = clean_sdf(sdf)
    _global.logger.write('info', "Converting to pdb using Openbabel")
    pdb = _ob(f'{_global.obabel} -isdf -opdb', stdin=sdf_text)
    
    return _common(sdf, pdb)
