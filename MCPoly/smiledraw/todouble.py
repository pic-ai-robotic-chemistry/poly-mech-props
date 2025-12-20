import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem
import os
import sys

def todouble(smiles, i, ringchoose='left', verbose=False):
    frags = split(smiles)
    lenfrags = len(frags)
    num = None
    num3 = None
    m = 0
    k = 0
    num0 = re.search(r'[0-9]+', frags[i-1])
    if num0:
        num = eval(num0.group(0))
    num0 = re.search(r'[0-9]+', frags[i])
    if num0:
        num3 = eval(num0.group(0))
    limits = ['O', 'S', 'F', 'Cl', 'Br', 'I']
    for limit in limits:
        if limit in frags[i]:
            raise TypeError("This element can't add double bonds.")
    if 'N' in frags[i] and '(' in frags[i+1]:
        raise TypeError("This structure can't add C=N now.")
    limits = ['c', 'n', 'o']
    for limit in limits:
        if limit in frags[i]:
            raise TypeError("Can't add double bond on conjugative system.")
    if '=' in frags[i] or '=' in frags[i-1]:
        raise TypeError('This part is already with a double bond.')
    if '#' in frags[i] or '#' in frags[i-1]:
        raise TypeError('This part is already with a triple bond.')
    if ')' in frags[i-1] and '(' in frags[i]:
        raise TypeError('Too much substitutions to add double bonds.')
    elif '(' in frags[i]:
        j = i
        while 1:
            if ')' in frags[j]:
                if '(' in frags[j+1]:
                    raise TypeError('Too much substitutions to add double bonds.')
                else:
                    break
            j = j + 1
    
    if '[C@@H]' in frags[i-1]:
        dt = re.split(r'\[C@@H\]', frags[i-1])
        frags[i-1] = ''
        frags[i-1] = dt[0] + 'C' + dt[1]
    elif '[C@@]' in frags[i-1]:
        raise TypeError('Too much substitutions to add double bonds.')
    elif '[C@H]' in frags[i-1]:
        dt = re.split(r'\[C@H\]', frags[i-1])
        frags[i-1] = dt[0] + 'C' + dt[1]
    elif '[C@]' in frags[i-1]:
        raise TypeError('Too much substitutions to add double bonds.')

    if '[C@@H]' in frags[i]:
        dt = re.split(r'\[C@@H\]', frags[i])
        frags[i] = ''
        frags[i] = dt[0] + 'C' + dt[1]
        if ')' in frags[i-1]:
            j = i - 1
            while 1:
                if '(' in frags[j]:
                    if ')' in frags[j-1]:
                        raise TypeError('Too much substitutions to add double bonds.')
                    elif '[C@@H]' in frags[j-1]:
                        dt = re.split(r'\[C@@H\]', frags[j-1])
                        frags[j-1] = dt[0] + 'C' + dt[1]
                    elif '[C@H]' in frags[j-1]:
                        dt = re.split(r'\[C@H\]', frags[j-1])
                        frags[j-1] = dt[0] + 'C' + dt[1]
                    else:
                        break
                j = j - 1
    elif '[C@@]' in frags[i]:
        raise TypeError('Too much substitutions to add double bonds.')
    elif '[C@H]' in frags[i]:
        dt = re.split(r'\[C@H\]', frags[i])
        frags[i] = ''
        frags[i] = dt[0] + 'C' + dt[1]
        if ')' in frags[i-1]:
            j = i - 1
            while 1:
                if '(' in frags[j]:
                    if ')' in frags[j-1]:
                        raise TypeError('Too much substitutions to add double bonds.')
                    elif '[C@@H]' in frags[j-1]:
                        dt = re.split(r'\[C@@H\]', frags[j-1])
                        frags[j-1] = dt[0] + 'C' + dt[1]
                    elif '[C@H]' in frags[j-1]:
                        dt = re.split(r'\[C@H\]', frags[j-1])
                        frags[j-1] = dt[0] + 'C' + dt[1]
                    else:
                        break
                j = j - 1
    elif '[C@]' in frags[i]:
        raise TypeError('Too much substitutions to add double bonds.')
    elif ')' in frags[i-1]:
        j = i - 1
        while 1:
            if '(' in frags[j]:
                if ')' in frags[j-1]:
                    raise TypeError('Too much substitutions to add double bonds.')
                elif '[C@@H]' in frags[j-1]:
                    dt = re.split(r'\[C@@H\]', frags[j-1])
                    frags[j-1] = dt[0] + 'C' + dt[1]
                elif '[C@H]' in frags[j-1]:
                    dt = re.split(r'\[C@H\]', frags[j-1])
                    frags[j-1] = dt[0] + 'C' + dt[1]
                else:
                    break
            j = j - 1
    
    if '(' in frags[i]:
        z1 = re.findall(r'\(+', frags[i])
    else:
        z1 = []

    if ringchoose == 'right' and num3 != None:
        num2 = re.split(r'[0-9]+', frags[i])
        print(frags[i])
        frags[i] = num2[0] + '=' + str(num3) + num2[1]
        if verbose == True:
            print(frags[i])
            print(frags)
    else:
        frags[i] = '(' * len(z1) + '=' + frags[i][len(z1):]
        if verbose == True:
            print(frags[i])
            print(frags)
    SMILES = merge(frags)
    mole_SMILES = Chem.MolFromSmiles(SMILES)
    #display(mole_SMILES)
    return Chem.MolToSmiles(mole_SMILES)

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

def merge(frags):
    x = ''
    for frag in frags:
        x = x + frag
    return x