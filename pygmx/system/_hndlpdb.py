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
    
def readLine(line):
    vals = {}
    vals['record'] = line[:6].strip()
    
    if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
        vals['content'] = line[6:]
        return vals
    
    vals['serial'] = line[6:11].strip()
    vals['name'] = line[12:16].strip()
    vals['altLoc'] = line[16]
    vals['resName'] = line[17:20].strip()
    vals['chainID'] = line[21]
    vals['resSeq'] = line[22:26].strip()
    vals['iCode'] = line[26]
    vals['x'] = line[30:38].strip()
    vals['y'] = line[38:46].strip()
    vals['z'] = line[46:54].strip()
    vals['occupancy'] = line[54:60].strip()
    vals['tempFactor'] = line[60:66].strip()
    vals['element'] = line[76:78].strip()
    vals['charge'] = line[78:80].strip()
    
    return vals
            
def editLine(line, **kwargs):
    vals = readLine(line)
    vals.update(kwargs)
    return writeLine(**vals)

def matchLine(line, **kwargs):
    vals = readLine(line)
    flg = False
    for key, val in kwargs.items():
        if vals[key] == val: flg = True
        else: flg = False
    
    return flg