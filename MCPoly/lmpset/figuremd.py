import re
import os
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

def multicreatemol(filename, turn, width, height):
    try:
        ind = open('{0}_dump.xyz'.format(filename), 'r+')
    except:
        raise FileNotFoundError('Perhaps your _dump.xyz file has some errors. '
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
    #print(main)
    raw_mol = Chem.MolFromXYZBlock(main)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    v = draw_with_spheres(conn_mol, width, height)
    return v

def figuremd(file, loc='./', num=0, width=300, height=300):
    """
    A method to see the current geometry structure and optimization trajectory of the ORCA Molecular Dynamics, powered by py3Dmol and ipywidgets package.
    
    figuremd(file, loc='./', num=0, width=300, height=300)
    num: The step of your convergence. The default is the origin structure.
    width, height: The size of your 3D geometry molecule strcuture. Default: 300x300.
    After forming the 3D geometry molecule strcuture, you can scroll to see other structures of relevant molecules.
    """
    figure = step(loc, file) - 1
    def turn(num):
        path = os.getcwd()
        os.chdir(loc)
        multicreatemol(file, num, width, height)
        os.chdir(path)
    interact(turn, num=iw.IntSlider(min=0, max=figure, step=0, value=num))