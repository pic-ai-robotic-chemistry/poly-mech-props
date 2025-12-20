import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem

def Sring_selection(smiles, connections, verbose=False):
    #e.g. connections = {2:'R', 8:'S'} 
    frags = split(smiles)
    lenfrags = len(frags)
    i = list(connections.keys())
    types = list(connections.values())
    i1 = i[0]
    i2 = i[1]
    num = 0
    for j, frag in enumerate(frags):
        num0 = re.search(r'\+?\-?[0-9]+', frag)
        if num0:
            num1 = re.search(r'\+', frag)
            if num1:
                continue
            num2 = re.search(r'\-', frag)
            if num2:
                continue
            if eval(num0.group(0)) > num:
                num = eval(num0.group(0))
    num = num + 1
    if verbose == True:
        print(num)
    for j, frag in enumerate(frags):
        if j == i1:
            status = 0
        elif j == i2:
            status = 1
        if j == i1 or j == i2:
            if types[status] != '1' and 'C' not in frag:
                raise TypeError("This is not the carbon atom, it can only be set '1'.")
            right_bracket = re.search(r'\)+', frag)
            if right_bracket:
                right_bracket_num = len(right_bracket.group(0))
                frags[j] = frag[:-right_bracket_num] + str(num) + ')' * right_bracket_num
            else:
                frags[j] = frag + str(num)
            xx = re.search(r'\[C\@+\]', frags[j])
            if xx:
                raise TypeError("These two atoms can't be connected.")
            parts_a = re.split(r'\[C\@H\]',frags[j])
            if verbose == True:
                print(parts_a)
            if len(parts_a) >= 2:
                if types[status] == '1':
                    frags[j] = parts_a[0] + '[C@]' + parts_a[1]
                    continue
                else:
                    raise TypeError("The pattern of connection can only be set '1'.")
            parts_a = re.split(r'\[C\@\@H\]',frags[j])
            if verbose == True:
                print(parts_a)
            if len(parts_a) >= 2:
                if types[status] == '1':
                    frags[j] = parts_a[0] + '[C@@]' + parts_a[1]
                    continue
                else:
                    raise TypeError("The pattern of connection can only be set '1'.")
            parts_a = re.split(r'[A-Z]',frags[j])
            if verbose == True:
                print(1, parts_a, frags[j])
            branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j, verbose=False)
            if branch[0] != '' and branch[1] != '':
                raise TypeError("These two atoms can't be connected.")
            if len(parts_a) >= 2:
                if types[status] == '1':
                    pass
                elif types[status] == 'R':
                    if branch[0] == '':
                        frags[j] = parts_a[0] + '[C@H]' + parts_a[1]
                    elif branch[1] == '':
                        frags[j] = parts_a[0] + '[C@]' + parts_a[1]
                elif types[status] == 'S':
                    if branch[0] == '':
                        frags[j] = parts_a[0] + '[C@@H]' + parts_a[1]
                    elif branch[1] == '':
                        frags[j] = parts_a[0] + '[C@@]' + parts_a[1]
            elif len(parts_a) == 1:
                if types[status] == '1':
                    pass
                elif types[status] == 'R':
                    if branch[0] == '':
                        frags[j] = '[C@H]' + parts_a[0]
                    elif branch[1] == '':
                        frags[j] = '[C@]' + parts_a[0]
                elif types[status] == 'S':
                    if branch[0] == '':
                        frags[j] = '[C@@H]' + parts_a[0]
                    elif branch[1] == '':
                        frags[j] = '[C@@]' + parts_a[0]
            
    if verbose == True:
        print(frags)
    SMILES = merge(frags)
    mole_SMILES = Chem.MolFromSmiles(SMILES)
    #xx = Chem.MolToSmiles(mole_SMILES)
    #parts = re.split('[CH]', xx)
    #result = ''
    #for i,part in enumerate(parts):
    #    if i == 0:
    #        result = result + part
    #    else:
    #        result = 'C' + result + part
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

def bracket_module(frags, j, verbose=False):
    branch = ['', '']
    brack_left = 0
    brack_right = 0
    branch_frag1 = 0
    branch_frag2 = 0
    k = 0
    l = 0
    while 1:
        try:
            if '(' in frags[j+1]:
                z1 = re.findall(r'\(',frags[j+1])
                z2 = re.findall(r'\)',frags[j+1])
                brack_left = brack_left + len(z1)
                brack_right = brack_right + len(z2)
                branch[k] = branch[k] + frags[j+1]
                j = j + 1
                l = 1
                if k == 0:
                    branch_frag1 = branch_frag1 + 1
                elif k == 1:
                    branch_frag2 = branch_frag2 + 1
            elif l == 1:
                z1 = re.findall(r'\(',frags[j+1])
                z2 = re.findall(r'\)',frags[j+1])
                brack_left = brack_left + len(z1)
                brack_right = brack_right + len(z2)
                branch[k] = branch[k] + frags[j+1]
                j = j + 1
                if k == 0:
                    branch_frag1 = branch_frag1 + 1
                elif k == 1:
                    branch_frag2 = branch_frag2 + 1
            else:
                break
            if brack_left == brack_right:
                l = 0
                k = k + 1
                if '(' in frags[j+1]:
                    pass
                else:
                    break
        except:
            break
    if verbose == True:
        print(branch)
    return branch, brack_left, brack_right, branch_frag1, branch_frag2