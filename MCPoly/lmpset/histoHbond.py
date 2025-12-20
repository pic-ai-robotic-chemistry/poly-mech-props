import os
import sys
from ase.io import read
from ase import Atoms
import numpy as np

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from Hbond import Hbond

mydir = os.path.dirname( __file__ )
lmpdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(lmpdir)
from figuremdsave import figuremdsave

def step(loc, file, fmt):
    path = os.getcwd()
    os.chdir(loc)
    if fmt == 'dump':
        file = file + '_dump.xyz'
    elif fmt == 'xyz':
        file = file + '.xyz'
    f = open(file, 'r')
    for i,line in enumerate(f):
        if i == 0:
            atomnum = eval(line[:-1])
    f.close()
    os.chdir(path)
    return (i + 1) // (atomnum + 2)

def histoHbond(file, loc='./', between='', HBond=2.0, angle=150, mode='all', fmt='dump'):
    path = os.getcwd()
    os.chdir(loc)
    figure = step(loc, file, fmt)
    if between == '':
        between = range(figure)
    AllHBonds = []
    for i in between:
        figuremdsave(file, loc='./', num=i, save='intermediate', fmt=fmt)
        HBonds = Hbond('intermediate', loc, HBond=HBond, angle=angle, mode=mode)
        AllHBonds.append(len(HBonds))
        os.system('rm intermediate.xyz')
    os.chdir(path)
    return AllHBonds