# This is part of MiMiCPy

"""

This module contains the Script class that serves
as the base class for all other script. It takes care
of assining parameters with the dot operator.

"""

from collections import OrderedDict 
from ..utils.errors import ScriptError
from .._global import _Global as gbl
from abc import ABC, abstractmethod

class Script(ABC):
    
    def __init__(self, **kwargs):
        # create an internal orderdict to store all params and retain order
        super().__setattr__('__orddict__', OrderedDict())
        
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def __setattr__(self, key, value):
        if key.startswith('_') or key == 'hasparams' or key == 'params':
            # protected/private attrs and hasparams() and params() are stored in the normal __dict__
            self.__dict__[key] = value
        else:
            # rest are script params and stored in __orddict__
            # all params are lower case
            self.__orddict__[key.replace(' ', '--').replace('-', '_').lower()] = value
        
    def __getattr__(self, key):
        if key.startswith('_') or key == 'hasparams' or key == 'params':
            # if asking fro protected/private attrs, hasparams() or params() look in the normal __dict__
            return self.__getattribute__('__dict__')[key]
        else:
            key = key.lower()
            # else look in __orddict__
            if key not in self.__getattribute__('__orddict__'):
                raise ScriptError(key) # raise param not found error
                
            return self.__getattribute__('__orddict__')[key]
    
    def hasparam(self, key):
        """Convienience function to check if script has a param"""
        if key in self.__orddict__:
            return True
        else:
            return False
    
    def params(self):
        """Convienience function to return dict of params"""
        return self.__orddict__
    
    def __repr__(self):
        """For printty printing in Jupyter"""
        return self.__str__()
    
    @classmethod    
    def fromFile(cls, script):
        return cls.fromText(gbl.host.read(script))
                
    @classmethod
    @abstractmethod  
    def fromText(cls, text):
        # should be defined in children
        pass
    
    @abstractmethod  
    def __str__(self):
        # should be defined in children
        pass