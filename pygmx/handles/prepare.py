# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

#from mimicpy.parsers import pdb as hpdb
from .base import BaseHandle
from .._global import _Global as _global

class Prepare(BaseHandle):
    """
    Prepare MM topology, by running pdb2gmx, editconf, solvate, etc.
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, status=[]):
        """Class constructor"""
        
        # parameters to pass to gmx genion
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''}
        self._box_kwargs = {'c': '', 'd': 1.9, 'bt': 'cubic'} # parameters to pass to gmx editconf
        self._solavte_kwargs = {'cs': 'spc216.gro'} # parameters to pass to gmx solvate
        self._topol_kwargs = {'water': 'tip3p', 'ff': 'amber99sb-ildn'} # parameters to pass to gmx pdb2gmx
        self.his_str = '' # string version of list of histidine protonation states, input to pdb2gmx
        super().__init__(None,status) # call BaseHandle constructor to init _status dict
        
        if self._status['prep'] == '':
            self.dir = 'prepare' # dir of handle, can be changed by user
        else:
            self.dir = self._status['prep']
            
        _global.logger.write('debug', f'Set handle directory as {self.dir}..')
        
        # all files names used when running gmx, can be changed by user
        self.confin = "confin.pdb"
        self.conf = 'conf.pdb'
        self.conf1 = 'conf1.gro'
        self.conf2 = 'conf2.gro'
        self.conf3 = 'conf3.gro'
        self.topol = "topol.top"
        self.mpt = "topol.mpt"
        self.preproc = "topol.pptop"
        self.ions = "ions"
    
    def pdb2gmxParams(self, **kwargs):
        """
        
        Set self._topol_kwargs, custom parameters for pdb2gmx
        Also, stringifies histiidine protonation states list
        
        """
        
        # check if his = ['D1', 'E2', ...] was passed in kwargs
        if 'his' in kwargs:
            self.his_str = '\n'.join(kwargs['his']) # join list as string
            # the list contains protonation states as human readable D1, E2
            # convert to input gromacs expects, D1-> 0; E2-> 1
            self.his_str = self.his_str.replace('D1', '0')
            self.his_str = self.his_str.replace('E2', '1')
            kwargs['his'] = '' # set the his option on, so -his will be passed to pdb2gmx
        
        self._topol_kwargs = kwargs # set custom topol kwargs
    
    def editconfParams(self, **kwargs): self._box_kwargs = kwargs
        
    def prepProtein(self, protein):
        """Runs gmx pdb2gmx for pure protien, and adds ligands/waters to pdb and .top in the right order"""
        
        self.setcurrent(key='prepMM') # set the current prepMM directory to self.dir in _status dict
        
        _global.logger.write('info', 'Preparing protein topology..')
        
        _global.logger.write('debug2', f"Writing {protein.name} to {self.confin}..")
        
        _global.host.write(protein.pdb+protein.water, f'{self.dir}/{self.confin}')
        
        _global.logger.write('info', 'Calculating protein topology')
        
        # call pdb2gmx, depending on if his was set
        if self.his_str == '':
            _global.logger.write('debug', "Letting Gromacs calculate histidine protonantion states..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs)
        else:
            _global.logger.write('debug', "Reading histidine protonation states from list..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs, stdin=self.his_str)
        
        conf = f'{self.dir}/{self.conf}'
        
        pdb = _global.host.read(conf) # output of pdb2gmx
        
        _global.logger.write('debug2',  f"Output of pdb2gmx saved in {self.conf}")
        
        ######
        # Followinng section reads output of pdb2gmx, which contains only protein+water residues
        # extracts crystal structure water residues
        # combine protine residues+ligand+water (IN THAT ORDER!!) are rewrite output of pdb2gmx
        ######
        lines = []
        splt = pdb.splitlines()
        for i, line in enumerate(splt[::-1]):
            vals = hpdb.readLine(line) # parse PDB using system._hndlpdb functions
            if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
                lines.append(line) 
            elif vals['record'] == 'HETATM' or vals['record'] == 'ATOM':
                if vals['resName'] == 'HOH': lines.append(line)
        
        _global.logger.write('info', "Combining ligand structure and topology with protein..")
        # combine protein, ligand and water in that order, so that gromacs doesn't complain about order
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + protein.ligand_pdb + '\n'.join(lines[::-1])
        
        _global.host.write(conf_pdb, conf)
        
        self.gmx('editconf', f = self.conf, o = self.conf1, dirc=self.dir, **self._box_kwargs)
        
        ######
        # Followinng section combines .itp data of all ligands into a single file
        # gromacs doesn't likes many #include for many ligands
        # so we have to combine all into a single ligands.itp file
        # Note: we need to combine all [ atomtypes ] section of itp
        # and write it as the first section of ligand.itp
        # protien.ligands is an OrderedDict so order in which we add ligands to topology
        # (i.e., order in which we loop), is also the order in which the ligands were added to pdb above
        # this is imp!! the orders have to be the same
        ######
        
        # init atomtypes section of itp
        top1 = (f"; Topology data for all non-standard resiudes in {protein.name} created by MiMiCPy\n"
            "; AmberTools was used to generate topolgy parameter for Amber Force Field, conversion to GMX done using Acpype"
                    "\n\n[ atomtypes ]\n")
        top2 = '' # init rest of itp
        
        _global.logger.write('debug2', "Combining ligands topology into single file ligands.itp..")
        for ligname, lig in protein.ligands.items():
            t1, t2 = lig.splitItp() # splits ligand itp into contents of [ atomtypes ], and eveything else
	
            top1 += t1 # update atomtypes section
            top2 += t2 # update rest of itp
            
            _global.logger.write('debug2', f"Writing position restraint file for {lig.name}")
            
            _global.host.write(lig.posre, f"{self.dir}/posre_{lig.name}.itp") # write position restraint file
            
            # add ligand to [ molecule ] section, assumed to be last section of .top
            _global.host.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        _global.host.write(top1+'\n'+top2, f"{self.dir}/ligands.itp")
        
        topol = f"{self.dir}/{self.topol}"
        
        # add #include "ligands.itp" after #include "..../forcefield.itp"
        _global.host.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" '+topol)
        # SOL is in between protein and ligand in [ molecule ]
        # we need to remove that and add it to the end, to match pdb order
        _global.host.run('grep -v SOL topol.top > topol_.top && mv topol_.top '+topol)
        _global.host.run(f'echo SOL {protein.hoh_mols} >> '+topol)
        _global.logger.write('debug', "ligands.itp added to topol.top")
        
        _global.logger.write('info', 'Topology prepared..')
        
        nonstd_atm_types = {}
        for name, lig in protein.ligands.items():
            nonstd_atm_types.update( dict(zip(lig.atm_types, lig.elems)) )
        
        #generate MPT file
        self.getMPT(ligs=list(protein.ligands.values()))
        
        _global.logger.write('info', "Converting to gro..")
        
        self.gmx('editconf', f = self.conf, o = self.conf1, dirc=self.dir, **self._box_kwargs)
        
        self.saveToYaml()
        
    def genionParams(self, **kwargs): self._ion_kwargs = kwargs
    def solvateParams(self, **kwargs): self._solavte_kwargs = kwargs
    
    def getBox(self, genion_mdp):
        """Run gmx solvate, gmx grompp and gmx genion in that order"""
        
        _global.logger.write('info', "Solvating box..")
        self.gmx('solvate', cp = self.conf1, o = self.conf2, p = self.topol, dirc=self.dir, **self._solavte_kwargs)
        
        _global.logger.write('info', "Adding ions to neutralize charge..")
        self.grompp(genion_mdp, self.ions, gro = self.conf2, pp = self.preproc, dirc=self.dir)
        
        # sent SOL to stdin, so gromacs replaces some SOL molecules with ions 
        self.gmx('genion', s = self.ions_tpr, o = self.conf3, p = self.topol, dirc=self.dir, **self._ion_kwargs, stdin="SOL")
        
        _global.logger.write('info', 'Simulation box prepared..')
        
        self.saveToYaml()
    
    def prepSystem(self, sys_dir=None, ligs=None, mpt=None, guess_elems=False):
        """Get the MPT topology, used in prepare.QM"""
        
        if sys_dir.strip() != '' or sys_dir.strip() != None:
            self.dir = sys_dir
        
        if ligs==None: ligs = []
        
        # get non std residue dict from ligs
        nonstd_atm_types = {}
        for lig in ligs:
            nonstd_atm_types.update( dict(zip(lig.atm_types, lig.elems)) )
        
        # generate proprocessed topology for now, TO DO: changed writer mpt to read from .top
        self.grompp(mdp.MDP.defaultGenion(), self.ions, gro = self.conf2, pp = self.preproc, dirc=self.dir) 
        
        if mpt == None: mpt = f"{self.dir}/{self.mpt}" # if no mpt file was passed, use default value
        else: mpt = f"{self.dir}/{mpt}"
        
        mptwrite(self.preproc, mpt, nonstd_atm_types, guess=guess_elems)
        
        self.toYaml()
