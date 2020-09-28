# This is part of MiMiCPy

"""

This module contains the Protein class
to handle proteins

"""

from . import ligand as lig, _hndlWater
from .._global import _Global as gbl
from collections import OrderedDict

class Protein:
   
   def __init__(self, pdb, name):
        self.name = name
        
        gbl.logger.write('debug', f"Creating protein {name}..")
        
        self.pdb = pdb
        
        gbl.logger.write('debug', "Extracting water..")
        
        self.water = gbl.host.run('grep ^HETATM', stdin=self.pdb)
        self.water = gbl.host.run('grep HOH', stdin=self.water)
        
        gbl.logger.write('debug', "Extracting amino acid residues..")
        
        self.pdb  = f"TITLE     {self.name}\n" + gbl.host.run('grep ^ATOM', stdin=self.pdb)
        
        self.ligands=OrderedDict()
        
        self.ligand_pdb = ""
        
        self._lig_elems = [] # for MiMiC
        
   def addLigand(self, ligand):
        
       if not isinstance(ligand, lig.NonStdLigand):
           raise TypeError(f"Cannot add {type(ligand)} as ligand")
       
       if not isinstance(ligand, lig.StdLigand):
           gbl.logger.write('debug', f"Adding non standard ligand {ligand.name} to {self.name}..") 
           
           self.ligands.update({ligand.name:ligand})
           
           self.ligand_pdb += ligand.pdb
           
           self._lig_elems.extend(ligand.elems) # for MiMiC
           
       else:
           gbl.logger.write('debug', f"Adding standard ligand {ligand.name} to {self.name}..") 
           
           self.pdb += ligand.pdb
    
   def addLigands(self, ligand_list):
       for ligand in ligand_list:  self.addLigand(ligand)
       
   def stripWater(self, *args):
       # args: list of ['lig_name', 'chain', dist] eg. ['AKG', 'A', 3], ['AKG', 'B', 3]
      ids = [ r for (res, chain, dist) in args\
                             for r in _hndlWater.coordsWithin(self.ligands[res].pdb, chain, self.water, dist) ]
     
      self.hoh_mols = len(ids)
      
      command = '\|'.join(ids)
      
      self.water = gbl.host.run(f"grep '{command}'", stdin=self.water)
       
   @classmethod
   def fromFile(cls, pdb):
        gbl.host.checkFile(pdb)
        f = gbl.host.read(pdb)
        return cls(f, pdb.split('.')[0])
