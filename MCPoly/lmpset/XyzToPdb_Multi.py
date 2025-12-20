import re
import warnings

def XyzToPdb_Multi(file, size, steps=[], pattern=[90, 90, 90], loc='./', saveloc='./', savename=''):
    """
    The method to convert .xyz trajectory file into .pdb single file and trajectory file.
    XyzToPdb_Multi(file, size=[X, Y, Z], steps=[], pattern=[90, 90, 90], loc='./', saveloc='./', savename='')
    file: Your .pdb file name.
    size: The size of the box, including the x, y, z size.
    pattern: The structure of the cell.
    steps: The step number
    loc: File Location. The default is your current location.
    saveloc: Input File Save Location. The default is your current location.
    savename: The name of your saved .xyz file. Default is the current name of .pdb file.
    Example 1:
        Input:
            from MCPoly.lmpset import XyzToPdb_Multi
            XyzToPdb_Multi('Moletraj1', steps=[20], size=[100, 100, 100])
        Output in Moletraj1.pdb:
            AUTHOR    MCPoly
            REMARK    Step 1
            CRYST1  100.000  100.000  100.000  90.00  90.00  90.00
            ATOM      1 H                   52.027  49.421  55.325  0.00  0.00           H
            ATOM      2 O                   54.200  47.007  55.661  0.00  0.00           O
            ATOM      3 C                   51.875  48.349  55.484  0.00  0.00           C
            ......
            END
            REMARK    Step 2
            CRYST1  100.000  100.000  100.000  90.00  90.00  90.00
            ATOM      1 H                   52.016  49.419  55.328  0.00  0.00           H
            ATOM      2 O                   54.188  47.005  55.669  0.00  0.00           O
            ATOM      3 C                   51.864  48.348  55.487  0.00  0.00           C
            ......
            END
            With all structure in PDB file.
            
    Example 2:
        Input:
            from MCPoly.lmpset import XyzToPdb_Multi
            XyzToPdb_Multi('Moletraj1', steps=[20], size=[100, 100, 100])
        Output in Moletraj1.pdb:
            AUTHOR    MCPoly
            REMARK    Step 20
            CRYST1  100.000  100.000  100.000  90.00  90.00  90.00
            ATOM      1 H                   52.009  49.420  55.337  0.00  0.00           H
            ATOM      2 O                   54.188  47.012  55.673  0.00  0.00           O
            ATOM      3 C                   51.860  48.348  55.493  0.00  0.00           C
            ......
            

    Example 3:
        Input:
            from MCPoly.lmpset import XyzToPdb_Multi
            XyzToPdb_Multi('Moletraj1', steps=[-2, -1], size=[100, 100, 100])
        Output in Moletraj1:
            REMARK    Step 24
            CRYST1  100.000  100.000  100.000  90.00  90.00  90.00
            ATOM      1 H                   52.005  49.421  55.340  0.00  0.00           H
            ATOM      2 O                   54.186  47.013  55.674  0.00  0.00           O
            ATOM      3 C                   51.857  48.348  55.495  0.00  0.00           C
            ......
            END
            REMARK    Step 25
            CRYST1  100.000  100.000  100.000  90.00  90.00  90.00
            ATOM      1 H                   52.008  49.419  55.338  0.00  0.00           H
            ATOM      2 O                   54.186  47.009  55.674  0.00  0.00           O
            ATOM      3 C                   51.859  48.346  55.493  0.00  0.00           C
            ......
            Trajectory of the last two steps.
    """
    step = 1
    f = open(loc+file+'.xyz', 'r')
    if savename == '':
        savename = file
    try:
        w = open(saveloc+savename+'.pdb', 'x')
    except:
        w = open(saveloc+savename+'.pdb', 'w')
    stepnum = 1
    i = 0
    xmin = 99999
    xmax = -99999
    ymin = 99999
    ymax = -99999
    zmin = 99999
    zmax = -99999
    w.write('AUTHOR    MCPoly\n')
    j = 0
    for line in f:
        if j == 0:
            num = eval(line[:-1])
        if j in [0, 1]:
            j = j + 1
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
        elif len(b) == 0:
            stepnum = stepnum + 1
            j = 0
        j = j + 1
    j = 0
    f.close()
    f = open(loc+file+'.xyz', 'r')
    for line in f:
        if j == 1:
            if steps == [] or (step in steps) or (step-stepnum-1 in steps):
                w.write('REMARK    Step {0}\n'.format(step))
                w.write('CRYST1  {0:>7.3f}  {1:>7.3f}  {2:>7.3f}  {3:.2f}  {4:.2f}  {5:.2f}\n'.format(*size, *pattern))
        if j in [0, 1]:
            j = j + 1
            continue
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
            if steps == [] or (step in steps) or (step-stepnum-1 in steps):
                w.write('ATOM  {0:>5} {1:<2}                 '.format(j-1, ele.group(0)))
                w.write('{0:>7.3f} {1:>7.3f} {2:>7.3f} '.format(eval(b[0])+size[0]/2-(xmax+xmin)/2, eval(b[1])+size[1]/2-(ymax+ymin)/2, eval(b[2])+size[2]/2-(zmax+zmin)/2))
                w.write('{0:>5.2f} {1:>5.2f}           {2}\n'.format(0.00, 0.00, ele.group(0)))
                if j == num + 1:
                    w.write('END\n')
        elif len(b) == 0:
            step = step + 1
            j = 0
        j = j + 1
    f.close()
    i = 0
    w.close()