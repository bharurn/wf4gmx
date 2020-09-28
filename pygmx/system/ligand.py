# This is part of MiMiCPy

"""

This module contains NonStdLigand and StdLigand classes to load ligands and protnate/get topology

"""

from . import _hndlpdb as hpdb
from . import _addH, _getItp
from .._global import _Global as _global
from ..utils.errors import EnvNotSetError

class NonStdLigand: 
    
    def __init__(self, pdb, itp, posre, chains, name, elems, atm_types):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
        self.elems = elems # for MPT
        self.atm_types = atm_types # for MPT
    
    def splitItp(self):
        start = False
        val = ''
        val2 = ''
        for line in self.itp.splitlines():

            if line.strip() == '[ atomtypes ]':
                start = True
            elif start:
			    # stop reading atomtypes section when blank line or moleculetype is reached
                if line.strip() == '' or  line.strip() == '[ moleculetype ]':
                    start = False
                else: val += line + '\n'
		
            elif not start: # if not start, read line into val2
                val2 += line + '\n'
		
        return val, val2	

    def _matchpdb2itp(self):
        start = False
        vals = {}
        i = 1
        
        print("Matching atom order in pdb to itp..")

        print("Reading itp atoms..")
        for line in self.itp.splitlines(): 
            if "[ atoms ]" in line:
                start = True
            elif start:
                splt = line.split()
                if len(splt) == 0:
                    break
                if splt[0].isdigit():
                    vals[splt[4]] = i
                    i = i + 1
        
        pdb_str = ""
        natms = 0
        
        print("Renumbering pdb atoms..")
        for line in self.pdb.splitlines():
            splt = hpdb.readLine(line)
            if splt['record'] == "HETATM" or splt['record'] == "ATOM":
                no = vals[splt['name']]
                
                if no > natms and splt['chainID'] == 'A': natms = no
        
                no_str = str( no + natms*( ord(splt['chainID'])-ord('A') ) )
        
                line = hpdb.editLine(line, serial=no_str)
        
            pdb_str += line+'\n'
        
        self.pdb = pdb_str.replace('ATOM  ','HETATM')
        print("Sorting pdb atoms..")
        self.pdb = _global.host.run('sort -nk2', stdin=self.pdb)
    
    @classmethod
    def _load(cls, mol, pH, prep2pdb):
        pdb, chains, mol_name, elems = _addH.do(mol, pH)
        
        itp, posre, atm_types = _getItp.do(mol_name, prep2pdb)
        
        lig = cls(pdb, itp, posre, chains, mol_name, elems, atm_types)
        
        lig._matchpdb2itp()
        
        return lig
    
    @classmethod    
    def fromBlock(cls, mol, pH=7, prep2pdb={}):
        print(f"Creating non standard ligand from SDF data..")
        
        return cls._load(mol, pH, prep2pdb)
    
    @classmethod    
    def fromFile(cls, file, pH=7, prep2pdb={}):
        _global.host.checkFile(f'{file}.sdf')
        
        mol = _global.host.read(f'{file}.sdf')
        
        print(f"Creating non standard ligand from {file}..")
        
        return cls._load(mol, pH, prep2pdb)
    
    @staticmethod    
    def gaussFromSDF(mol, nc, pH=7, parallel=False):
        
        print(f"Running Gaussian calculation on ligand..")
        
        print(f"Converting to PDB..")
        
        pdb, chains, name, elems = _addH.do(mol, pH)
        
        _global.host.write(pdb, f"{name}.pdb")
        
        print("Generating Gaussian input file using AmberTools Antechamber..")
        
        _global.host.run(f"antechamber -i {name}.pdb -fi pdb -o {name}.com -fo gcrt -nc {nc}")
        
        _global.host.query_rate = 30
        _global.host.redirectStdout(f"{name}.out") # add this to local _global.host
        
        if _global.gauss is None:
            raise EnvNotSetError("No _global.gaussian executable given!")
        
        _global.host.runbg(f"{_global.gauss} {name}.com {name}.out", query_rate=0)
        
        print("Gaussian run submitted as a background job..\n"
              "Do not close host and/or this script untill the run is complete!!\n"
              f"Please check {name}.out for status..")
    
    @staticmethod
    def getPrepGauss(mol):
        print(f"Converting _global.gaussian output file to Amber prep..")
        _global.host.checkFile(f"{mol}.out")
        print('Running AmberTools antechamber..')
        out = _global.host.run(f"{_global.antechamber} -i {mol}.out -fi gout -o {mol}.prep -fo prepi -c resp -rn {mol}")
        print(f"Antechamber output dumped...\nResidue {mol} created in {mol}.prep..")
        
        return out
        
class StdLigand(NonStdLigand):
    @classmethod
    def fromBlock(cls, mol):
        print(f"Creating standard ligand from SDF block..")
        
        pdb, chains, mol_name, elems = _addH.donoH(mol)
        
        return cls(pdb, "", "", chains, mol_name, elems, {})
    
    @classmethod 
    def fromFile(cls, file):
        _global.host.checkFile(f'{file}.sdf')
        
        mol = _global.host.read(f'{file}.sdf')
        
        print(f"Creating standard ligand from {file}..")
        
        pdb, chains, mol_name, elems = _addH.donoH(mol)
        
        return cls(pdb, "", "", chains, mol_name, elems, {})
    
