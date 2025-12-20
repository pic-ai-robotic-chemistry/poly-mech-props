import re
import warnings

def XyzToPdb(file, size, pattern=[90, 90, 90], loc='./', saveloc='./', savename=''):
    """
    The method to convert .xyz file into .pdb file.
    XyzToPdb(file, size=[X, Y, Z], pattern=[90, 90, 90], loc='./', saveloc='./', savename='')
    file: Your .pdb file name.
    size: The size of the box, including the x, y, z size.
    pattern: The structure of the cell.
    loc: File Location. The default is your current location.
    saveloc: Input File Save Location. The default is your current location.
    savename: The name of your saved .xyz file. Default is the current name of .pdb file.
    TIPS: If there are some atoms out of the box, the code will warn.
    Example:
        Input:
            from MCPoly.lmpset import XyzToPdb
            XyzToPdb('Mole1', size=[100, 100, 100])
        Output in Mole1.xyz:
            AUTHOR    MCPoly
            CRYST1   100.000   100.000   100.000  90.00  90.00  90.00
            ATOM      1 H                   48.430  44.708  53.282  0.00  0.00           H
            ATOM      2 O                   47.934  45.536  53.413  0.00  0.00           O
            ATOM      3 C                   47.140  45.369  54.569  0.00  0.00           C
            ......
    """
    step = 1
    f = open(loc+file+'.xyz', 'r')
    if savename == '':
        savename = file
    try:
        w = open(saveloc+savename+'.pdb', 'x')
    except:
        w = open(saveloc+savename+'.pdb', 'w')
    stepnum = 0
    i = 0
    xmin = 99999
    xmax = -99999
    ymin = 99999
    ymax = -99999
    zmin = 99999
    zmax = -99999
    w.write('AUTHOR    MCPoly\n')
    w.write('CRYST1  {0:>7.3f}  {1:>7.3f}  {2:>7.3f}  {3:.2f}  {4:.2f}  {5:.2f}\n'.format(*size, *pattern))
    for j,line in enumerate(f):
        if j in [0, 1]:
            continue
        b = re.findall(r'\-?[0-9]+\.[0-9]+', line)
        if len(b) == 3:
            if eval(b[0]) < xmin:
                xmin = eval(b[0])
            if eval(b[0]) > xmax:
                xmax = eval(b[0])
            if eval(b[1]) < ymin:
                ymin = eval(b[1])
            if eval(b[1]) > ymax:
                ymax = eval(b[1])
            if eval(b[2]) < zmin:
                zmin = eval(b[2])
            if eval(b[2]) > zmax:
                zmax = eval(b[2])
    f.close()
    f = open(loc+file+'.xyz', 'r')
    for j,line in enumerate(f):
        if j in [0, 1]:
            continue
        i = i + 1
        ele = re.search(r'[A-Z][a-z]*', line)
        b = re.findall(r'\-?[0-9]+\.[0-9]+', line)
        #print(ele,b)
        if ele and len(b) == 3:
            if eval(b[0])+size[0]/2-(xmax+xmin)/2 < 0.0:
                warnings.warn("Some atoms are out of the cell box.")
            elif eval(b[1])+size[1]/2-(ymax+ymin)/2 < 0.0:
                warnings.warn("Some atoms are out of the cell box.")
            elif eval(b[2])+size[2]/2-(zmax+zmin)/2 < 0.0:
                warnings.warn("Some atoms are out of the cell box.")
            elif eval(b[0])+size[0]/2-(xmax+xmin)/2 > size[0]:
                warnings.warn("Some atoms are out of the cell box.")
            elif eval(b[1])+size[1]/2-(ymax+ymin)/2 > size[1]:
                warnings.warn("Some atoms are out of the cell box.")
            elif eval(b[2])+size[2]/2-(zmax+zmin)/2 > size[2]:
                warnings.warn("Some atoms are out of the cell box.")
            w.write('ATOM  {0:>5} {1:<2}                 '.format(i, ele.group(0)))
            w.write('{0:>7.3f} {1:>7.3f} {2:>7.3f} '.format(eval(b[0])+size[0]/2-(xmax+xmin)/2, eval(b[1])+size[1]/2-(ymax+ymin)/2, eval(b[2])+size[2]/2-(zmax+zmin)/2))
            w.write('{0:>5.2f} {1:>5.2f}           {2}\n'.format(0.00, 0.00, ele.group(0)))
            #print('{0:>7.3f} {1:>7.3f} {2:>7.3f} '.format(eval(b[0])+size[0]/2+(xmax+xmin)/2, eval(b[1])+size[1]/2+(ymax+ymin)/2, eval(b[2])+size[2]/2+(zmax+zmin)/2))
    f.close()
    w.write('END\n')
    i = 0
    w.close()