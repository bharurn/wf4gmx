# This is part of MiMiCPy

"""

This module contains the Slurm class that
allows for pythonic creation/manipulation of
Slurm jobscripts

"""

from .base import Script

class Slurm(Script):
    def __init__(self, name='jobscript', shebang='/bin/bash', cmd_hdr='srun', cmds = [], **kwargs):
        #if kwargs is None:
        #    kwargs = {'nodes':2, 'ntasks': 16, 'cpus_per_task': 3, 'mem_per_cpus': 3, 'export': 'NONE', 'time': '00:10:00'}
            
        super().__init__(**kwargs)
        
        self._shebang = shebang
        self._name = name
        self._cmds = cmds
        self._cmd_hdr = cmd_hdr
    
    def setSpecial(self, name=None, shebang=None, cmd_hdr=None):
        if name: self._name = name
        if shebang: self._shebang = shebang
        if cmd_hdr: self._cmd_hdr = cmd_hdr
        
    @classmethod    
    def fromText(cls, text):
        shebang = "/usr/sh"
        name = 'jobscript'
        cmds = []
        kwargs = {}
        for line in text.splitlines():
            if line.startswith('#!'):
                  shebang = line[2:].strip()
            elif line.startswith('#SBATCH -A'):
                  kwargs.update({'account':line[11:].strip()})
            elif line.startswith('#SBATCH --'):
                l = line[10:].replace('-', '_').split('=')
                key = l[0].strip()
                val = l[1].strip()
                if key == 'job_name':
                    name = val
                elif key == 'output':
                    continue
                elif key != 'name':
                    kwargs.update({key:val})
            elif line.strip() != '':
                cmds.append(line)
        
        return cls(name=name, shebang=shebang, cmds=cmds, **kwargs)

    def add(self, cmd, **kwargs):
        opt = ' '
        for k,v in kwargs.items():
            if len(k) == 1: opt += f"-{k}{v} "
            else: opt += f"--{k.replace('_','-')}={v} "
        self._cmds.append(f"{self._cmd_hdr}{opt}{cmd}")
    
    def addMany(self, cmds):
        self._cmds.extend(cmds)
    
    def __getattr__(self, val):
        if val == 'job_name' or val == 'name':
            return self._name
        elif val == 'output':
            return f'{self._name}.%J.out'
        else:
            return super().__getattr__(val)
    
    def noCommands(self):
        if len(self._cmds) <= 0:
            return True
        return False
    
    def clearCommands(self, keep_source=False, keep_module=False):
        if keep_source:
            s = []
            for i, c in enumerate(self._cmds):
                if c.split()[0] == 'source':
                    s.append(self._cmds[i])
            self._cmds = s.copy() 
        elif keep_module:
            s = []
            for i, c in enumerate(self._cmds):
                if c.split()[0] == 'module':
                    s.append(self._cmds[i])
            self._cmds = s.copy()
        else:            
            self._cmds = []
        
    def __str__(self):
        cmd = f'#!{self._shebang}\n'
        
        for d in self.params():
            if getattr(self,d) == None: continue
            d_ = d.replace('_','-')
            cmd += f"#SBATCH --{d_}={getattr(self, d)}\n"
        
        cmd += f'#SBATCH --job-name={self._name}\n'
        cmd += f'#SBATCH --output={self._name}.%J.out\n\n'
        
        cmd += '\n'.join(self._cmds)
        
        return cmd