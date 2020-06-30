# This is part of MiMiCPy

"""

This module contains handles to run MD and MiMiC simulations

"""


from .base import BaseHandle
from .._global import _Global as _global
from ..utils.errors import EnvNotSetError
from ..utils.constants import hartree_to_ps
from . import _qmhelper

class MD(BaseHandle):
    """
    Simulte MD runs using grompp, convetr-tpr, mdrun
    Inherits from .core.base.BaseHandle
    
    """
    def mdrun(self, new, **kwargs):
        """Execute gmx mdrun"""
        
        def _do(new, **kwargs):
            """Runs mdrun, depending on if jobscript is there or not"""
            if self.jobscript is None:
                self.gmx('mdrun', **kwargs, deffnm = new)
            else:
                self.jobscript.add(self.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True))
        
        if 'noappend' in kwargs: # set noappend if noappend=True passed
            if kwargs['noappend'] == True:
                del kwargs['noappend']
                _do(new, **kwargs, noappend = '')
            else: _do(new, **kwargs)
        else:
            _do(new, **kwargs)
        
        if self.jobscript:
            if 'dirc' in kwargs: # for jobscripts, dirc has to be passed to sbatch() not gmx()
                dirc = kwargs['dirc']
                del kwargs['dirc']
            else:
                dirc = ''
                    
            jid = _global.host.sbatch(self.jobscript, dirc=dirc)
            _global.logger.write('info', "Gromacs simulations submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            _global.logger.write('info', "Gromacs simulation is now running as a background job..")
            if not _global.host.isLocal():
                _global.logger.write('info', "Please do not close remote host until the job is done!")
    
    def restart(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        """
        Restart a simulation, either:
            - from crash, in which case mdrun is called
            - increase the timesteps of a completed run by extending run time until time t in ps,
              in this case convert-tpr has to called with -until
            - increase timesteps of a completed run by adding time t in ps to current time,
            in this case convert-tpr has to called with -extend
        """
        
        self.setcurrent(new)
        
        if fromcrash:
            out = self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
        elif until != 0:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
        elif extend != 0:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = '{cpt}.cpt', noappend = noappend, dirc=new)
        
        self.saveToYaml()
        return out
        
    def run(self, mdp, **kwargs):
        """
        Wrapper around gromp() and mdrun()
        Also takes care of directory creation
        """
        new = mdp.name.lower().replace(' ', '_')
        
        if 'dir' in kwargs: # dir is an argument for the directory, if dir = None, then use cwd
            if kwargs['dir'] == None: kwargs['dir'] = ''
            _dir = kwargs['dir']
        else:   
            _dir = new
        self.setcurrent(_dir)
        
        _global.logger.write('info', f"Starting classical MD calculation: {new}..")
        
        _global.logger.write('debug', f"All files will be saved to the directory {_dir}..")
        
        self.grompp(mdp, f'{new}', dirc=_dir, **kwargs)
        
        out = self.mdrun(f'{new}', dirc=_dir)
        
        self.saveToYaml()
        
        return out
        
class MiMiC(BaseHandle):
    """
    Runs MiMiC runs by running both gmx mdrun and cpmd
    Inherits from .core.base.BaseHandle
    
    """
    
    def __cpmd(self, inp, out, onlycmd=False, dirc=''):
        """
        Function to execute cpmd command in "pythonic" way
        e.g., to execute cpmd cpmd.in path/to/pp > cpmd.oput
        call cpmd('cpmd.in', 'cpmd.out')
        Make sure cpmd_pp is set in _Global before that!
        """
        
        # check for env variables
        if _global.cpmd is None or _global.cpmd.strip() == '':
            raise EnvNotSetError('CPMD executable', 'cpmd')
        
        if _global.cpmd_pp is None or _global.cpmd_pp.strip() == '':
            raise EnvNotSetError('CPMD pseudopotential', 'cpmd_pp')
        
        cmd = f"{_global.cpmd} {inp} {_global.cpmd_pp} > {out}"
        
        if onlycmd: return cmd # return only cmd
        
        _global.logger.write('debug', cmd)
        
        _global.logger.write('debug', "Running {cmd}..")
        _global.host.runbg(cmd, dirc=dirc)
    
    # functions to set individual slurm parameters for gmx and cpmd
    def setGMXSettings(self, **kwargs): self.gmx_opt = kwargs
    def setCPMDSettings(self, **kwargs): self.cpmd_opt = kwargs
    
    def run(self, mdp, inp, tpr=None, dirc=''):
        """
        Writes inp to CPMD file
        And runs gmx mdrun and cpmd depending
        Also supports Slurm
        """
        new = inp.info.lower().replace(' ', '_')
        self.setcurrent(new)
        _global.host.mkdir(f"{new}/cpmd") # cpmd
        _global.host.mkdir(f"{new}/gmx") # and gmx directories
        
        ######create mimic tpr file
        _global.logger.write('debug2', "Changing Gromacs integrator to MiMiC..")
        mdp.integrator = 'mimic'
    
        _global.logger.write('debug', f"Writing atoms in QM region to {self.index}..")
        QMMM_grps = 'QMatoms'
        
        # set no of steps in mdp file
        if not mdp.hasparam('nsteps'): mdp.nsteps = 1000 # default value, if not present in mdp
        # set timestep in mdp file
        if not mdp.hasparam('dt'): mdp.dt = 0.001 # default value, if not present in mdp
        
        # write index file
        _global.host.write(_qmhelper.index(inp._ndx, self.QMMM_grps), _global.host.join(self.dir, self.index))
        mdp.QMMM_grps = QMMM_grps
        
        _global.logger.write('info', "Generating Gromacs TPR file for MiMiC run..")
        
        self.grompp(mdp, 'mimic.tpr', gro=self.gro, n=self.index, dirc=f'{new}/gmx')
        ######
        
        ####write cpmd file
        inp.mimic.paths = f"1\n{_global.host.pwd()+new}/gmx" # set path in cpmd script
        # set no of steps
        inp.cpmd.maxsteps = mdp.nsteps
        # set timestep
        inp.cpmd.timestep = round(mdp.dt/hartree_to_ps) # convert to atomic units and round off
        
        _global.logger.write('debug2', "Set PATH in MIMIC section as {inp.mimic.paths}")
        _global.host.write(str(inp), f"{new}/cpmd/{new}.inp")
        ####
        
        # below is similar procedure to simulate.MD.run()
        if self.jobscript:
            if not hasattr(self, 'gmx_opt'): self.gmx_opt = {}
            if not hasattr(self, 'cpmd_opt'): self.cpmd_opt = {}
            
            self.jobscript.add(self.gmx('mdrun', deffnm=f'gmx/mimic', onlycmd=True, dirc='new'), **self.gmx_opt)
            self.jobscript.add(self.__cpmd(f"cpmd/{new}.inp", f"cpmd/{new}.out", onlycmd=True, dirc='new'), **self.cpmd_opt)
            jid = _global.host.sbatch(self.jobscript, dirc=dirc)
            _global.logger.write('info', "MiMiC run submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            self.gmx('mdrun', deffnm=f'gmx/mimic', dirc='new')
            self.__cpmd(f"cpmd/{new}.inp", f"cpmd/{new}.out", dirc='new')
            _global.logger.write('info', "MiMiC simulation is now running in the background..")
            
            if not _global.host.isLocal():
                _global.logger.write('info', "Please do not close remote host until the job is done!")
            
        self.saveToYaml()
        
        