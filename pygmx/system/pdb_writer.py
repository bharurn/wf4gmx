# This is part of MiMiCPy

"""

This module contains the helper functions to efficiently
handle pdb lines by converting them to a dictionary

"""

from mimicpy.parser.pdb import readLine

def _getSpaces(s, n):
    return (n - len(s))*' '

def _errHandle(s, mx):
    s = s.strip()
    if len(s) > mx:
        return s[:len(s)-mx+1]
    elif len(s) < 1:
        return ' '
    else:
        return s

def writeLine(record='ATOM', serial='1', name='XX', altLoc=' ', resName='XXX', chainID='A',\
              resSeq='1', iCode=' ', x='00.00', y='00.00', z='00.00', occupancy='1.00', tempFactor='00.00',
              element='C', charge=' '):
    
    record = _errHandle(record, 6)
    record = record + _getSpaces(record, 6)
    
    serial = _errHandle(serial, 5)
    serial = _getSpaces(serial, 5) + serial
    
    name = _errHandle(name, 4)
    #name = ' ' + name
    name = _getSpaces(name, 4) + name
    
    altLoc = _errHandle(altLoc, 1)
    
    resName = _errHandle(resName, 3)
    resName = resName + _getSpaces(resName, 3)
    
    chainID = _errHandle(chainID, 1)
    
    resSeq = _errHandle(resSeq, 1)
    resSeq = _getSpaces(resSeq, 4) + resSeq
    
    iCode = _errHandle(iCode, 1)
    
    x = _errHandle(x, 8)
    x = _getSpaces(x, 8) + x
    y = _errHandle(y, 8)
    y = _getSpaces(y, 8) + y
    z = _errHandle(z, 8)
    z = _getSpaces(z, 8) + z
    
    occupancy = _errHandle(occupancy, 6)
    occupancy = _getSpaces(occupancy, 6) + occupancy
    
    tempFactor = _errHandle(tempFactor, 6)
    tempFactor = _getSpaces(tempFactor, 6) + tempFactor
    
    element = ''.join([e for e in element if e.isalpha()])
    element = _errHandle(element, 2)
    element = _getSpaces(element, 2) + element
    
    charge = _errHandle(charge, 2)
    charge = _getSpaces(charge, 2) + charge
    
    return record + serial + ' ' + name + altLoc + resName + ' ' + chainID + resSeq + iCode + ' '*3 +\
            x + y + z + occupancy + tempFactor + ' '*10 + element + charge

def editLine(line, **kwargs):
    vals = readLine(line)
    vals.update(kwargs)
    return writeLine(**vals)