import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem

def S_C_selection(substitution, smiles, i, verbose=False):
    '''This function will change the carbon into other substitution.'''
    frags = split(smiles)
    lenfrags = len(frags)
    Cs = ['O', 'N', 'S', 'C(=O)', 'C(=N)', 'S(=O)', 'S(=O)(=O)']
    Cs_special = ['O', 'S', 'C(=O)', 'C(=N)', 'S(=O)', 'S(=O)(=O)']
    Cs_double = ['N', 'n']
    keyfrag = frags[i]
    l0 = [i,8]
    l = l0[0] + 1
    o = 0
    x = 0

    double_check = 0
    m = l0[0]
    n = 0
    if '=' in keyfrag:
        double_check = 2
    else:
        try:
            m = m + 1
            if '(' in frags[m]:
                n = 1
            while ')' in frags[m] or n == 1:
                m = m + 1
                n = 0
            if '=' in frags[m]:
                double_check = 1
        except:
            pass

    q = 1
    if '(' in frags[i][0] and '(' in frags[i+1][0]:
        raise TypeError("This carbon atom can't be replaced.")
    if '(' in frags[i][0] and '(' in frags[i+1][0]:
        raise TypeError("This carbon atom can't be replaced.")
    if '=' in frags[i+1][0] and substitution in Cs_special:
        raise TypeError("This site can only be replace by N.")
    if '=' in frags[i][0] and substitution in Cs_special:
        raise TypeError("This site can only be replace by N.")
    if '(' in frags[i+q][0]:
        while 1:
            bracket_left = re.findall(r'\(', frags[i+q])
            o = o + len(bracket_left)
            bracket_right = re.findall(r'\)', frags[i+q])
            o = o - len(bracket_right)
            if verbose == True:
                print(o)
            try:
                frags[i+q] == frags[i+q]
            except:
                break
            if o == 0:
                if '(' == frags[i+q+1][0]:
                    raise TypeError("This carbon atom can't be replaced.")
                elif '=' in frags[i+q+1]:
                    raise TypeError("This carbon atom can't be replaced.")
                else:
                    break
            q = q + 1
        if substitution in Cs_special:
            raise TypeError("This site can only be replace by N.")
    
    if 'C' in keyfrag and double_check != 0 and substitution in Cs_double:
        return C_double(substitution, frags, i, double_check, verbose=verbose)
    elif 'c' in keyfrag and substitution in Cs_double:
        return c_double(substitution, frags, i, verbose=verbose)
    elif 'C' in keyfrag and substitution in Cs:
        return C_normal(substitution, frags, i, verbose=verbose)
    else:
        raise TypeError("This substitution can't be set.")

def addfrag(frags, i, verbose=False):
    frags.append('')
    j = len(frags) - 2
    while 1:
        if j >= i:
            frags[j+1] = frags[j]
            if verbose == 'True':
                print(frags)
            j = j - 1
        else:
            break
    frags[i+1] = ''
    if verbose == True:
        print('Final:',frags)
    return frags

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

def C_normal(substitution, frags, i, verbose=False):
    z1 = re.match(r'\(', frags[i])
    z2 = re.search(r'\)', frags[i])
    num0 = re.search(r'[0-9]+', frags[i])
    if z1:
        left = z1.group(0)
    else:
        left = ''
    if z2:
        right = z2.group(0)
    else:
        right = ''
    if num0:
        num = num0.group(0)
    else:
        num = ''
    if left != '' and num != '':
        if substitution in ('S', 'O', 'S(=O)', 'S(=O)(=O)'):
            raise TypeError("This carbon atom can't" + 
                            "be replaced by {0}.".format(substitution))
    frags[i] = left + substitution + num + right
    SMILES = merge(frags)
    return SMILES

def C_double(substitution, frags, i, double_check, verbose=False):
    z1 = re.match(r'\(', frags[i])
    z2 = re.search(r'\)', frags[i])
    num0 = re.search(r'[0-9]+', frags[i])
    if z1:
        left = z1.group(0)
    else:
        left = ''
    if z2:
        right = z2.group(0)
    else:
        right = ''
    if num0:
        num = num0.group(0)
    else:
        num = ''
    if double_check == 1:
        frags[i] = left + substitution + num + right
    elif double_check == 2:
        frags[i] = '=' + substitution + num + right
    SMILES = merge(frags)
    return SMILES

def c_double(substitution, frags, i, verbose=False):
    z1 = re.match(r'\(', frags[i])
    z2 = re.search(r'\)', frags[i])
    num0 = re.search(r'[0-9]+', frags[i])
    if z1:
        raise TypeError("This carbon atom can't be replaced.")
    if z2:
        right = z2.group(0)
    else:
        right = ''
    if num0:
        num = num0.group(0)
    else:
        num = ''
    frags[i] = substitution + num + right
    SMILES = merge(frags)
    return SMILES