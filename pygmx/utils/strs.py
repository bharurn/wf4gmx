# This is part of MiMiCPy

"""

This module contains the f function to simulate
f-string behavior in python version 3.6-

"""

import inspect
import re

def f(s):
    frame = inspect.currentframe().f_back
    v = dict(**frame.f_globals)
    v.update(**frame.f_locals)
    return s.format(s, **v)

def clean(txt, comment=None):
    if comment:
        txt = re.sub(re.compile(f"{comment}(.*)\n" ) ,"\n" , txt) # strip comments
    return re.sub(re.compile("\n{2,}" ) ,"\n" , txt) # remove double new lines