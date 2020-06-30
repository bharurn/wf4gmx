# This is part of MiMiCPy

"""

This module contains classes for all custom exceptions
used in the package

"""

class MiMiCPyError(Exception):
   """Generic exception from MiMiCPy"""

class SelectionError(MiMiCPyError):
    """Error in selecting atoms in MPT"""
    
class ParserError(MiMiCPyError):
    """
    Raised when a file could not be parsed
    using the format requsted
    """
    def __init__(self, file='data', ftype='', no='', extra=''):
        self.file = file
        self.ftype = ftype
        self.no = no
        self.extra = extra
        
        if ftype: self.ftype = ' as '+ ftype
        if self.no: self.no = f" at line number {no}"
        if self.extra: self.extra = ': '+ extra
        
    def __str__(self):
        s = f"Error parsing {self.file}{self.ftype}{self.no}{self.extra}"
        
        return s

class ExecutionError(MiMiCPyError):
    """
    Raised when a command executed by host exits
    with and error message
    """
    def __init__(self, cmd, msg):
        self.cmd = cmd
        self.msg = msg
        
    def __str__(self):
        return f"Command attempted {self.cmd}\n{self.msg}"

class GromacsError(ExecutionError):
    """
    Raised when a Gromacs command executed by host
    exits with and error message
    """

class SlurmBatchError(ExecutionError):
    """
    Raised when a Slurm sbatch command executed by host
    exits with and error message
    """
    def __str__(self):
        return f"Attempted to run sbatch {self.cmd}\nself.msg"

class EnvNotSetError(MiMiCPyError):
    """
    Raised when an enviornment path requested by
    MiMiCPy is not set
    """
    def __init__(self, env, cmd):
        self.env = env
        self.cmd = cmd
        
    def __str__(self):
        return f"{self.env} not set! Please set it with the keyword {self.cmd}."

class ScriptError(MiMiCPyError):
    """
    Raised when a requested parameter of a script
    object does not exist
    """
    def __init__(self, msg):
        self.msg = msg
        
    def __str__(self):
        return f"{self.msg} not a paramater of the script"

def defaultHook(cmd, out):
    """Default error checking hook called by host"""
    if 'error' in out.lower() or 'not recognisable' in out.lower() or 'not found' in out.lower():
        raise ExecutionError(cmd.split(';')[-1], out)
    
def asserter(boolean, error, *args, **kwargs):
    """Custom assert to raise custom errors"""
    if not boolean: raise error(*args, **kwargs)