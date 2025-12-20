import os
import re
from ase import Atoms
from ase.io import write

def XyzToData(file, x=[-99999,99999], y=[-99999,99999], z=[-99999,99999], loc='./', savename=''):
    if savename == '':
        savename = file
    opath = os.getcwd()
    os.chdir(loc)
    elements = []
    elementsort = []
    elementdict = {}
    w1 = 0
    n1 = 0
    l = []
    f = open(file+'.xyz', 'r')
    try:
        w = open(loc+savename+'.data', 'x')
    except:
        w = open(loc+savename+'.data', 'w')
    w.write('#Powered by MCPoly\n\n')
    j = 1
    for i,line in enumerate(f):
        if i == 0:
            a1 = eval(line[:-1])
        if i >= 2:
            ele = re.search(r'[A-Z][a-z]?', line)
            ele = ele.group(0)
            elements.append(ele)
            if ele not in elementsort:
                elementsort.append(ele)
                elementdict[ele] = j
                j = j + 1
    f.close()
    w.write('{0} atoms\n'.format(a1))
    w.write('{0} atom types\n\n'.format(len(elementsort)))
    w.write('{0:>9.5f} {1:>9.5f} xlo xhi\n'.format(x[0], x[1]))
    w.write('{0:>9.5f} {1:>9.5f} ylo yhi\n'.format(y[0], y[1]))
    w.write('{0:>9.5f} {1:>9.5f} zlo zhi\n\n'.format(z[0], z[1]))
    w.write('Masses\n\n')
    for i,ele in enumerate(elementsort):
        if ele == 'H':
            w.write('{0} {1}\n'.format(i+1, 1.008))
        elif ele == 'C':
            w.write('{0} {1}\n'.format(i+1, 12.011))
        elif ele == 'N':
            w.write('{0} {1}\n'.format(i+1, 14.007))
        elif ele == 'O':
            w.write('{0} {1}\n'.format(i+1, 15.999))
        elif ele == 'S':
            w.write('{0} {1}\n'.format(i+1, 32.060))
        elif ele == 'F':
            w.write('{0} {1}\n'.format(i+1, 19.000))
        elif ele == 'Cl':
            w.write('{0} {1}\n'.format(i+1, 35.455))
    w.write('\nAtoms\n\n')
    j = 1
    f = open(file + '.xyz', 'r')
    for i,line in enumerate(f):
        if i >= 2:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            current_ele = elements[j-1]
            try:
                w.write('{0} {1} 0 {2:>9.5f} {3:>9.5f} {4:>9.5f}\n'.format(j, elementdict[current_ele],
                                                                          eval(l1[0]), eval(l1[1]), eval(l1[2])))
                j = j + 1
            except:
                pass
    w.write('\n')
    f.close()
    w.close()
    os.chdir(opath)
    print(elementdict)
    return elementdict