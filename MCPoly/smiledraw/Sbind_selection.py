import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem
import os
import sys

mydir = os.path.dirname( __file__ )
smilesdir = os.path.join(mydir, '..', 'smiledraw')
sys.path.append(smilesdir)
from smiledraw.Ssub_selection import Ssub_selection
from smiledraw.Sring_selection import Sring_selection

def Sbind_selection(connections, smiles, num, verbose=False):
    #e.g. connections = {2:'R', 8:'S'}
    if num < 1:
        raise TypeError("Error on binding number.")
    
    i = list(connections.keys())
    types = list(connections.values())
    a = Ssub_selection({types[0]: 'C' * num}, smiles, i[0], verbose=verbose, special=True)
    if verbose == True:
        print(a)
    b = Sring_selection(a, {i[0]+num: '1', i[1]+num: types[1]}, verbose=verbose)
    return b