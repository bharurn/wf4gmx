# This is part of MiMiCPy

"""

This module contains the Logger, LogString and StdOut
classes to allow for logging multiple streams and redirection

"""
from .._global import _Global as gbl
from .errors import MiMiCPyError

class Logger:
    def __init__(self, **kwargs):
       
        for k,v in kwargs.items():
            setattr(self, k, v)
    
    def add(self, **kwargs):
       
        for k,v in kwargs.items():
            setattr(self, k, v)
        
    def write(self, option, value):
        writer = getattr(self, option)
        if writer != None:
            if not hasattr(writer, 'write'):
                raise MiMiCPyError(f'{option} is not a logger stream')
            writer.write(value+'\n')
        
class LogString():
    def __init__(self):
        self.val = ''
    
    def write(self, s):
        self.val += s
    
    def __str__(self):
        return self.val
    
    def __repr__(self):
        return self.val
    
    def __getattr__(self, attr):
        # for all other string methods
        return getattr(self.val, attr)
    
    def __eq__(self, other):
        if isinstance(other, LogString):
            return self.val == other.val
        elif isinstance(other, str):
            return self.val == other
        else:
            return False
    
class LogFile:
    def __init__(self, name, forceLocal=False):
        self.fname = name
        self.forceLocal = forceLocal
        
    def write(self, s):
        if not self.forceLocal:
            gbl.host.write(s, self.fname)
        else:
            with open(self.fname, 'r') as f:
                f.write(s)