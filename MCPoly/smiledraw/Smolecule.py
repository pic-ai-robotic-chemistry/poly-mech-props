import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem
import sys
import os
import py3Dmol

mydir = os.path.dirname( __file__ )
smilesdir = os.path.join(mydir, '..', 'smiledraw')
sys.path.append(smilesdir)
from smiledraw.Ssub_selection import Ssub_selection
from smiledraw.Sbind_selection import Sbind_selection
from smiledraw.Sring_selection import Sring_selection
from smiledraw.S_C_selection import S_C_selection
from smiledraw.todouble import todouble

#import Ssub_selection
#import Sbind_selection
#import Sring_selection
#import S_C_selection

def smi2conf(mol):
    '''Convert SMILES to rdkit.Mol with 3D coordinates'''
    if mol is not None:
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol)
        return mol
    else:
        return None

def MolTo3DView(mol, size=(400, 550), surface=False, opacity=0.5):
    """Draw molecule in 3D (Not Origin)
    
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'cartoon'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    mblock = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({'sphere':{'radius':0.4}, 'stick':{'radius':0.1}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    Chem.MolToXYZFile(mol, 'Mole_Untitled.xyz')
    return viewer

def split(smiles):
    frags = []
    a = re.findall(r'\(*\=?\#?\[?\\?\/?[A-Z]?[a-z]?\+?\-?[0-9]*\@*H?\]?[0-9]*\)*', smiles)
    for i, va in enumerate(a):
        if [a[i-2], a[i-1], a[i]] == ['([N+]', '([O-]', '=O)']:
            del frags[-1]
            del frags[-1]
            frags.append('([N+]([O-])=O)')
        elif [a[i-2], a[i-1], a[i]] == ['[N+]', '([O-]', '=O']:
            del frags[-1]
            del frags[-1]
            frags.append('[N+]([O-])=O')
        elif [a[i-2], a[i-1], a[i]] == ['[N+]', '([O-]', '=O)']:
            del frags[-1]
            del frags[-1]
            frags.append('[N+]([O-])=O)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(S', '(=O)', '(=O)', 'O)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('(S(=O)(=O)O)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['S', '(=O)', '(=O)', 'O']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('S(=O)(=O)O')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['S', '(=O)', '(=O)', 'O)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('S(=O)(=O)O)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(P', '(=O)', '(O)', 'O)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('(P(=O)(O)O)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['P', '(=O)', '(O)', 'O']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('P(=O)(O)O')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['P', '(=O)', '(O)', 'O)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('P(=O)(O)O)')
        elif [a[i-1], a[i]] == ['(C', '(=O))']:
            del frags[-1]
            frags.append('(C(=O))')
        elif [a[i-1], a[i]] == ['(C', '(=O)']:
            del frags[-1]
            frags.append('(C(=O)')
        elif [a[i-1], a[i]] == ['C', '(=O)']:
            del frags[-1]
            frags.append('C(=O)')
        elif [a[i-1], a[i]] == ['C', '(=O))']:
            del frags[-1]
            frags.append('C(=O))')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['C', '(F)', '(F)', 'F']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('C(F)(F)F')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(C', '(F)', '(F)', 'F)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('C(F)(F)F)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(C', '(F)', '(F)', 'F))']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('(C(F)(F)F)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['C', '(Cl)', '(Cl)', 'Cl']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('C(Cl)(Cl)Cl')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(C', '(Cl)', '(Cl)', 'Cl)']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('C(Cl)(Cl)Cl)')
        elif [a[i-3], a[i-2], a[i-1], a[i]] == ['(C', '(Cl)', '(Cl)', 'Cl))']:
            del frags[-1]
            del frags[-1]
            del frags[-1]
            frags.append('(C(Cl)(Cl)Cl)')
        elif a[i] == 'Cc':
            frags = addfrag(frags, i)
            frags.append('C')
            frags.append('c')
        elif a[i] == '(Cc':
            frags = addfrag(frags, i)
            frags.append('(C')
            frags.append('c')
        elif a[i] == '':
            continue
        else:
            frags.append(a[i])
    return frags

class Smolecule:
    def __init__(self, smiles):
        self.smiles = smiles
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)

    def S_sub(self, substitutions, i, verbose=False):
        self.smiles = Ssub_selection(substitutions, self.smiles, i, verbose=verbose)
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)
        return self.smiles

    def S_C(self, substitutions, i, verbose=False):
        self.smiles = S_C_selection(substitutions, self.smiles, i, verbose=verbose)
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)
        return self.smiles

    def S_bind(self, substitutions, i, verbose=False):
        self.smiles = Sbind_selection(substitutions, self.smiles, i, verbose=verbose)
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)
        return self.smiles

    def S_ring(self, connections, verbose=False):
        self.smiles = Sring_selection(self.smiles, connections, verbose=verbose)
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)
        return self.smiles

    def S_todouble(self, i, ringchoose='left', verbose=False):
        self.smiles = todouble(self.smiles, i, ringchoose=ringchoose, verbose=verbose)
        self.image = Chem.MolFromSmiles(self.smiles)
        self.frags = split(self.smiles)
        return self.smiles

    def view2d(self):
        return display(self.image)

    def view3d(self, size=(400,550)):
        conf = smi2conf(self.image)
        return MolTo3DView(conf,size).show()