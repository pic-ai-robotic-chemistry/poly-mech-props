import re
import os
import sys
import py3Dmol
from ipywidgets import interact
import ipywidgets as iw
from rdkit import Chem
from IPython import display
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDetermineBonds

def draw_with_spheres(mol, width, height):
    v = py3Dmol.view(width=width, height=height)
    IPythonConsole.addMolToView(mol, v)
    v.zoomTo()
    v.setStyle({'sphere':{'radius':0.4}, 'stick':{'radius':0.1}});
    v.show()
    return v

def setxyz(filename, savename, turn, fmt='dump'):
    try:
        if fmt == 'dump':
            ind = open('{0}_dump.xyz'.format(filename), 'r+')
        elif fmt == 'xyz':
            ind = open('{0}.xyz'.format(filename), 'r+')
    except:
        raise FileNotFoundError('Perhaps your _trj.xyz file has some errors. '
                                + 'Please check the location of '
                                + 'your XYZ file and then modify your file.')
    mains = ind.readlines()
    num = eval(mains[0][:-1])
    term = len(mains)/(num+2)
    main = ''
    partmain = mains[turn*(num+2):(turn+1)*(num+2)]
    for i, line in enumerate(partmain):
        if i == 1:
            main = main + '\n'
        else:
            main = main + partmain[i]
    if savename == 'intermediate':
        w = open('{0}.xyz'.format(savename), 'w')
    else:
        w = open('{0}_{1}.xyz'.format(savename, turn), 'w')
    w.write(main)
    ind.close()
    w.close()

def step(loc, file, choose=0):
    path = os.getcwd()
    os.chdir(loc)
    file = file + '_dump.xyz'
    f = open(file, 'r')
    for i,line in enumerate(f):
        if i == 0:
            atomnum = eval(line[:-1])
    f.close()
    os.chdir(path)
    return (i + 1) // (atomnum + 2)

def figuremdsave(file, loc='./', num=0, save='', fmt='dump'):
    """
    A method to save the current geometry structure of the ORCA optimisation.
        
    figure(num=0, save='')
    num: The step of your convergence.
    save: The name of the single status XYZ file you want to save. The default is not saving the single status.
    """
    if save == '':
        save = file
    try:
        setxyz(file, save, num, fmt)
    except:
        raise ValueError("The index 'num' is out of range.")