import math as m
import numpy as np
import re
from ase import Atoms
from rdkit import Chem

def Ssub_selection(substitutions, smiles, i, verbose=False, special=False):
    '''This function will change the hydrogen into other substitution.'''
    keys = list(substitutions.keys())
    values = list(substitutions.values())
    frags = split(smiles)
    lenfrags = len(frags)
    subs = ['(F)', '(Cl)', '(Br)', '(I)', '(O', '([N+]([O-])=O)', '(S', '(N']
    subs_b = ['F', 'Cl', 'Br', 'I', 'O', '[N+]([O-])=O', 'S', 'N']
    keyfrag = frags[i]
    l0 = [i,8]
    l = l0[0] + 1
    o = 0
    x = 0
        
    if 'H' in values and ['H'] != values:
        for i,key in enumerate(keys):
            if substitutions[key] == 'H':
                del substitutions[key]

    specialC1 = re.search(r'\[?C\@*H?\]?[0-9]+', keyfrag)
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
        
    if specialC1:
        return CringSub(substitutions, frags, i, specialC1.group(0),\
                        verbose=verbose, special=special)
    elif 'C' in keyfrag and double_check == 1:
        return CdoubleSub(substitutions, frags, i, 1, verbose=verbose, special=special)
    elif 'c' in keyfrag and double_check == 1:
        return CdoubleSub(substitutions, frags, i, 1, verbose=verbose, special=special)
    elif 'C' in keyfrag and double_check == 2:
        return CdoubleSub(substitutions, frags, i, 2, verbose=verbose, special=special)
    elif 'c' in keyfrag and double_check == 2:
        return CdoubleSub(substitutions, frags, i, 2, verbose=verbose, special=special)
    elif 'C' in keyfrag: #('H', 'F', 'Cl', 'Br', 'I', 'O', 'N', 'NO2', 'C', 'C=C', 'C#C', 'C=N', 'CHO', 'CN', 'SO3')
        return CSub(substitutions, frags, i, verbose=verbose, special=special)
    elif 'c' in keyfrag:
        return CSub(substitutions, frags, i, verbose=verbose, special=special)
    elif '[N+]' in keyfrag:
        return XtraNSub(substitutions, frags, i, verbose=verbose, special=special)
    elif 'N' in keyfrag:
        return NSub(substitutions, frags, i, verbose=verbose, special=special)
    elif 'O' in keyfrag or 'S' in keyfrag:
        return OSub(substitutions, frags, i, verbose=verbose, special=special)
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

def doublebondstyle(smiles, style, double_bond_no=-1, j=-1, output='frags'):
    frags = split(smiles)
    if double_bond_no != -1:
        n = 0
        for i,frag in enumerate(frags):
            if '=' in frag:
                n = n + 1
            if n == double_bond_no:
                Nunn = i
                j = i - 1
                k = 0
                m = None
                while 1:
                    if ')' not in frags[i-1] and j == i - 1:
                        m = 0
                        for l,char in enumerate(frags[j]):
                            if m == 0:
                                if char == '/' or char == '\\':
                                    m = char
                                    break
                                else:
                                    m = l
                                    break
                        if type(m) == int:
                            try:
                                frags[j] = frags[j][:m] + '/' + frags[j][m:]
                                m == '/'
                                break
                            except:
                                frags[j] = '/' + frags[j]
                                m == '/'
                                break
                    if ')' in frags[j]:
                        k = 1
                    if '(' in frags[j]:
                        m = 0
                        for l,char in enumerate(frags[j-1]):
                            if m == 0:
                                if char == '/' or char == '\\':
                                    m = char
                                    j = j - 1
                                    break
                                else:
                                    m = l
                                    j = j - 1
                                    break
                        if type(m) == int:
                            try:
                                frags[j] = frags[j][:m] + '/' + frags[j][m:]
                                m = '/'
                                break
                            except:
                                frags[j] = '/' + frags[j][m:]
                                m = '/'
                                break
                break
        try:
            i = Nunn + 1
            if i == Nunn + 1:
                if '/' in frags[i] or '\\' in frags[i]:
                    for char in frags[i]:
                        if m == '/':
                            if style == 'Z':
                                char = '/'
                            elif style == 'E':
                                char = '\\'
                            break
                        elif m == '\\':
                            if style == 'Z':
                                char = '\\'
                            elif style == 'E':
                                char = '/'
                            break
                else:
                    for c,char in enumerate(frags[i]):
                        if frags[i][0] != '(':
                            c = -1
                            break
                        if frags[i][c] == '(' and frags[i][c+1] != '(':
                            break
                    if c == -1:
                        if m == '/':
                            if style == 'Z':
                                frags[i] = '/' + frags[i]
                            elif style == 'E':
                                frags[i] = '\\' + frags[i]
                        if m == '\\':
                            if style == 'Z':
                                frags[i] = '\\' + frags[i]
                            elif style == 'E':
                                frags[i] = '/' + frags[i]
                    else:
                        if m == '/':
                            if style == 'Z':
                                frags[i] = frags[i][:c+1] + '/' + frags[i][c+1:]
                            elif style == 'E':
                                frags[i] = frags[i][:c+1] + '\\' + frags[i][c+1:]
                        if m == '\\':
                            if style == 'Z':
                                frags[i] = frags[i][:c+1] + '\\' + frags[i][c+1:]
                            elif style == 'E':
                                frags[i] = frags[i][:c+1] + '/' + frags[i][c+1:]
        except:
            pass
            
    elif j != -1:
        bond_num = 0
        for i,frag in enumerate(frags):
            if '=' in frag:
                bond_num = bond_num + 1
            if i == j:
                if '=' not in frag:
                    raise TypeError('There is no C=C bond in this site.')
                else:
                    return doublebondstyle(smiles, style, double_bond_no=bond_num, output=output)
    if output == 'frags':
        return frags
    elif output == 'smiles':
        return merge(frags)

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

def reverse(symbol):
    if symbol == '/':
        return '\\'
    elif symbol == '\\':
        return '/'
    elif symbol == '@':
        return '@@'
    elif symbol == '@@':
        return '@'

def bracket_module2(frags, j, side, verbose=False):
    branch = ['', '']
    brack_left = 0
    brack_right = 0
    branch_frag1 = 0
    branch_frag2 = 0
    l = 0
    m = j
    n = 0
    m = j + 1
    while side == 1:
        try:
            m = m + 1
            if '(' in frag[m]:
                n = 1
            elif n == -1:
                break
            if ')' in frag[m] or n == 1:
                if l == 0:
                    brack_left = brack_left + 1
                elif l == 1:
                    brack_right = brack_right + 1
                branch[l] = branch[l] + frag[m]
                n = 0
                if l == 1:
                    break
            if '=' in frag[m]:
                n = -1
                l = l + 1
        except:
            break
    m = j
    while side == 2:
        try:
            m = m - 1
            if ')' in frag[m]:
                n = 1
            if '(' in frag[m] or n == 1:
                m = m - 1
                brack_left = brack_left + 1
                branch[0] = frag[m] + branch[0]
                n = 0
                break
        except:
            break
    m = j + 1
    while side == 2:
        try:
            m = m + 1
            if '(' in frag[m]:
                n = 1
            if ')' in frag[m] or n == 1:
                m = m + 1
                brack_right = brack_right + 1
                branch[1] = branch[1] + frag[m]
                n = 0
                break
        except:
            break
    if verbose == True:
        print(branch)
    return branch, brack_left, brack_right

def CSub(substitutions, frags, i, verbose=False, special=False):
    keys = list(substitutions.keys())
    big_num = 0
    if len(keys) == 2 and ({*keys} != {'R', 'S'}, {*keys} != {'1', '1'}):
        raise TypeError("For two substitutions, it can only use '2'.")
    elif ('R' in keys or 'S' in keys) and i == 0:
        raise TypeError("For the first atom, it can only use '1' or '2'.")
    if 'c' in frags[i]:
        if len(keys) == 2:
            raise TypeError("For two substitutions in conjugative systems, it can only use '1'.")
        elif '2' in keys:
            raise TypeError("In this conjugative systems, it can only use '1'.")
        elif 'R' in keys:
            raise TypeError("In this conjugative systems, it can only use '1'.")
            #substitutions['1'] = substitutions['R']
            #del substitutions['R']
        elif 'S' in keys:
            raise TypeError("In this conjugative systems, it can only use '1'.")
            #substitutions['1'] = substitutions['S']
            #del substitutions['S']
    keys = list(substitutions.keys())
    j = i
    branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j, verbose=verbose)
    if verbose == True:
        print(branch)
    if '@@' in frags[i]:
        property_branch = ['S', 'R']
    elif '@' in frags[i]:
        property_branch = ['R', 'S']
    else:
        property_branch = ['1', '1']

    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))
    
    for key in keys:
        if key == 'Z' or key == 'E':
            raise TypeError("This carbon atom is on the C-C bond.")
        keyname = substitutions[key]
        if substitutions[key] == 'H':
            substitutions[key] = ''
        if substitutions[key] == 'NO2':
            substitutions[key] = '[N+]([O-])=O'
        if substitutions[key] == 'CHO':
            substitutions[key] = 'C(=O)'
        if substitutions[key] == 'COOH':
            substitutions[key] = 'C(=O)O'
        if substitutions[key] == 'CN':
            substitutions[key] = 'C#N'
        if substitutions[key] == 'SO3H':
            substitutions[key] = 'S(=O)(=O)O'
        if substitutions[key] == 'PO3H2':
            substitutions[key] = 'P(=O)(O)O'
        if substitutions[key] == 'CF3':
            substitutions[key] = 'C(F)(F)F'
        if substitutions[key] == 'CCl3':
            substitutions[key] = 'C(Cl)(Cl)Cl'
        if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
            substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
        if substitutions[key] == 'ring5'  or substitutions[key] == 'cyclo5':
            substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
        if substitutions[key] == 'ring6' or substitutions[key] == 'cyclo6':
            substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
        if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
            substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
        if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
            substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)
        try:
            frags[i+1] == frags[i+1]
        except:
            if key in ('R', 'S'):
                substitutions['1'] = substitutions[key]
                del substitutions[key]
                key = '1'
        a2_result = ''
        a = re.match(r'\(+', frags[i])
        a2 = re.search(r'\)+', frags[i])
        if a2:
            a2_result = a2.group(0)
        if a2 and key in ('R', 'S'):
            substitutions['1'] = substitutions[key]
            del substitutions[key]
            key = '1'
        if key == '1':
            if substitutions[key] == '':
                raise TypeError("Substitution 'H' can be only used in 'R', 'S' and '2'.")
            if a:
                if 'C' in frags[i]:
                    frags[i] = a.group(0) + 'C'
                elif 'c' in frags[i]:
                    frags[i] = a.group(0) + 'c'
                    for i in range(branch_frag1):
                        del frags[i+1]
                if branch[1] != '':
                    for i in range(branch_frag2):
                        del frags[i+branch_frag1+1]
                a_result = a.group(0)
            else:
                if 'C' in frags[i]:
                    frags[i] = 'C'
                elif 'c' in frags[i]:
                    frags[i] = 'c'
                    for j in range(branch_frag1):
                        del frags[i+1]
                if branch[1] != '':
                    for j in range(branch_frag2):
                        del frags[i+branch_frag1+1]
                a_result = ''
            frags = addfrag(frags, i+branch_frag1, verbose=verbose)
            if substitutions[key] == '':
                pass
            elif i+1 == len(frags)-1:
                frags[i+branch_frag1+1] = substitutions[key]
                i = i + 1
            else:
                frags[i+branch_frag1+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
                i = i + 1
        elif key == 'R':
            if a:
                if branch_frag1 == 0:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = a.group(0) + '[C@H]'
                elif '' in branch and '[C@@H' in frags[i]:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = a.group(0) + '[C@H]'
                elif '' in branch and '[C@H' in frags[i]:
                    for i in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + 'C'
                elif '' in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = a.group(0) + '[C@]'
                elif '' not in branch and '[C@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = a.group(0) + '[C@]'
                elif '' not in branch and '[C@@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = a.group(0) + '[C@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                elif '' not in branch and '[C@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = a.group(0) + '[C@]'
                    for j in range(branch_frag1):
                        del frags[i+1]
                elif '' not in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = a.group(0) + '[C@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = '[C@H]'
                elif '' in branch and '[C@@H' in frags[i]:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = '[C@H]'
                elif '' in branch and '[C@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = 'C'
                elif '' in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@H]'
                    else:
                        frags[i] = '[C@]'
                elif '' not in branch and '[C@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = '[C@H]'
                    else:
                        frags[i] = '[C@]'
                elif '' not in branch and '[C@@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@H]'
                    else:
                        frags[i] = '[C@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                elif '' not in branch and '[C@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@H]'
                    else:
                        frags[i] = '[C@]'
                    for j in range(branch_frag1):
                        del frags[i+1]
                elif '' not in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@H]'
                    else:
                        frags[i] = '[C@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            elif i+1 == len(frags)-1:
                frags[i+1] = substitutions[key] + len(a2_result) * ')'
                i = i + 1
            else:
                frags[i+1] = '(' + substitutions[key] + ')' + len(a2_result) * ')'
                i = i + 1
        elif key == 'S':
            if a:
                if branch_frag1 == 0:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = a.group(0) + '[C@@H]'
                elif '' in branch and '[C@H' in frags[i]:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = a.group(0) + '[C@@H]'
                elif '' in branch and '[C@@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + 'C'

                elif '' in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@@H]'
                    else:
                        frags[i] = a.group(0) + '[C@@]'
                elif '' not in branch and '[C@@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@@H]'
                    else:
                        frags[i] = a.group(0) + '[C@@]'
                elif '' not in branch and '[C@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@@H]'
                    else:
                        frags[i] = a.group(0) + '[C@@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                elif '' not in branch and '[C@@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@@H]'
                    else:
                        frags[i] = a.group(0) + '[C@@]'
                    for j in range(branch_frag1):
                        del frags[i+1]
                elif '' not in branch and '@@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = a.group(0) + '[C@@H]'
                    else:
                        frags[i] = a.group(0) + '[C@@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = '[C@@H]'
                elif '' in branch and '[C@H' in frags[i]:
                    if substitutions[key] == '':
                        break
                    else:
                        frags[i] = '[C@@H]'
                elif '' in branch and '[C@@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = 'C'

                elif '' in branch and '@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@@H]'
                    else:
                        frags[i] = '[C@@]'
                elif '' not in branch and '[C@@H' in frags[i]:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    if substitutions[key] == '':
                        frags[i] = '[C@@H]'
                    else:
                        frags[i] = '[C@@]'
                elif '' not in branch and '[C@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@@H]'
                    else:
                        frags[i] = '[C@@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                elif '' not in branch and '[C@@]' in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@@H]'
                    else:
                        frags[i] = '[C@@]'
                    for j in range(branch_frag1):
                        del frags[i+1]
                elif '' not in branch and '@@' not in frags[i]:
                    if substitutions[key] == '':
                        frags[i] = '[C@@H]'
                    else:
                        frags[i] = '[C@@]'
                    for j in range(branch_frag2):
                        del frags[i+1+branch_frag1]
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            elif i+1 == len(frags)-1:
                frags[i+1] = substitutions[key] + len(a2_result) * ')'
                i = i + 1
            else:
                frags[i+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
                i = i + 1
        elif key == '2':
            if a:
                frags[i] = a.group(0) + 'C'
                a_result = a.group(0)
            else:
                frags[i] = 'C'
                a_result = ''
            for j in range(branch_frag1 + branch_frag2):
                del frags[i+1]
            frags = addfrag(frags, i, verbose=verbose)
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            else:
                frags[i+1] = '('+ substitutions[key] +')'
                if keyname == 'ring3' or keyname == 'cyclo3':
                    substitutions[key] = 'C{0}CC{0}'.format(big_num+2)
                if keyname == 'ring5' or keyname == 'cyclo5':
                    substitutions[key] = 'C{0}CCCC{0}'.format(big_num+2)
                if keyname == 'ring6' or keyname == 'cyclo6':
                    substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+2)
                if keyname == 'Benzene' or keyname == 'benzene' or keyname == 'Ph':
                    substitutions[key] = 'c{0}ccccc{0}'.format(big_num+2)
                if keyname == 'Naphthalene' or keyname == 'naphthalene' or keyname == 'NA':
                    substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+3, big_num+4)
                if i+2 == len(frags)-1:
                    frags[i+2] = substitutions[key] + len(a2_result) * ')'
                else:
                    frags[i+2] = '('+ substitutions[key] +')' + len(a2_result) * ')'
            
    SMILES = merge(frags)
    if verbose == True:
        print(SMILES)
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)

def CdoubleSub(substitutions, frags, i, side, verbose=False, special=False):
    key = list(substitutions.keys())[0]
    if verbose == True:
        print(side)
    if key == 'R' or key == 'S':
        raise TypeError("This carbon atom is on the C=C bond.")
    elif ('R' in keys or 'S' in keys) and i == 0:
        raise TypeError("For the first atom, it can only use '1' or '2'.")
    j = i
    branch, len_left, len_right = bracket_module2(frags, j, side)
    if verbose == True:
        print(branch)

    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))
    
    keyname = substitutions[key]
    if substitutions[key] == 'H':
        substitutions[key] = ''
    if substitutions[key] == 'NO2':
        substitutions[key] = '[N+]([O-])=O'
    if substitutions[key] == 'CHO':
        substitutions[key] = 'C(=O)'
    if substitutions[key] == 'COOH':
        substitutions[key] = 'C(=O)O'
    if substitutions[key] == 'CN':
        substitutions[key] = 'C#N'
    if substitutions[key] == 'SO3H':
        substitutions[key] = 'S(=O)(=O)O'
    if substitutions[key] == 'PO3H2':
        substitutions[key] = 'P(=O)(O)O'
    if substitutions[key] == 'CF3':
        substitutions[key] = 'C(F)(F)F'
    if substitutions[key] == 'CCl3':
        substitutions[key] = 'C(Cl)(Cl)Cl'
    if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
        substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
    if substitutions[key] == 'ring5' or substitutions[key] == 'cyclo5':
        substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
    if substitutions[key] == 'ring6'or substitutions[key] == 'cyclo6':
        substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
    if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
        substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
    if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
        substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)
    try:
        frags[i+1] == frags[i+1]
        if key == '2':
            raise TypeError("In this structure, it can only add one substitution.")
    except:
        pass
    if '' in branch and (key == 'Z' or key == 'E'):
        substitutions['1'] = substitutions[key]
        del substitutions[key]
        key = '1'
    a2_result = ''
    a = re.match(r'\(+', frags[i])
    a2 = re.search(r'\)+', frags[i])
    if a2:
        a2_result = a2.group(0)
    if verbose == True:
        print(substitutions)
    if key in ('1', 'Z', 'E'):
        if a:
            if 'C' in frags[i]:
                if side == 1:
                    frags[i] = a.group(0) + 'C'
            elif 'c' in frags[i]:
                if side == 1:
                    frags[i] = a.group(0) + 'c'
            if side == 1:
                for j in range(len_left):
                    del frags[i+1]
            elif side == 2:
                for j in range(len_right):
                    del frags[i+1]
            a_result = a.group(0)
        else:
            if 'C' in frags[i]:
                if side == 1:
                    frags[i] = 'C'
                elif side == 2:
                    frags[i] = '=C'
            elif 'c' in frags[i]:
                if side == 1:
                    frags[i] = 'c'
                elif side == 2:
                    frags[i] = '=c'
            if side == 1:
                for j in range(len_left):
                    del frags[i+1]
            elif side == 2:
                for j in range(len_right):
                    del frags[i+1]
            a_result = ''
        
        frags = addfrag(frags, i, verbose=verbose)

        if substitutions[key] == '':
            pass
        elif i+1 == len(frags)-1:
            frags[i+1] = substitutions[key]
        else:
            frags[i+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
    
    elif key == '2':
        proper_branch = ''
        if a:
            if 'C' in frags[j]:
                if side == 1:
                    frags[j] = a.group(0) + 'C'
            elif 'c' in frags[j]:
                if side == 1:
                    frags[j] = a.group(0) + 'c'
            if side == 1:
                ring_num = 0
                for j in range(len_left):
                    s = re.search(r'[A-Z][a-z]?[0-9]+', frags[i+j+1])
                    if s:
                        ring_num = ring_num + 1
                    del frags[i+1]
                if ring_num % 2 == 1:
                    raise TypeError("This structure is on the ring, "
                                    + "so it can't change two substitutions.")
            elif side == 2:
                j = 0
                ring_num = 0
                while 1:
                    try:
                        s = re.search(r'[A-Z][a-z]?[0-9]+', frags[i+j+1])
                        if s:
                            ring_num = ring_num + 1
                        del frags[i+1]
                    except:
                        break
                    j = j + 1
                if ring_num % 2 == 1:
                    raise TypeError("This structure is on the ring, "
                                    + "so it can't change two substitutions.")
            a_result = a.group(0)
        else:
            if 'C' in frags[i]:
                if side == 1:
                    frags[i] = 'C'
                elif side == 2:
                    frags[i] = '=C'
            elif 'c' in frags[i]:
                if side == 1:
                    frags[i] = 'c'
                elif side == 2:
                    frags[i] = '=c'
            if side == 1:
                ring_num = 0
                for j in range(len_left):
                    s = re.search(r'[A-Z][a-z]?[0-9]+', frags[i+j+1])
                    if s:
                        ring_num = ring_num + 1
                    del frags[i+1]
                if ring_num % 2 == 1:
                    raise TypeError("This structure is on the ring, "
                                    + "so it can't change two substitutions.")
            elif side == 2:
                j = 0
                ring_num = 0
                while 1:
                    try:
                        s = re.search(r'[A-Z][a-z]?[0-9]+', frags[i+j+1])
                        if s:
                            ring_num = ring_num + 1
                        del frags[i+1]
                    except:
                        break
                    j = j + 1
                if ring_num % 2 == 1:
                    raise TypeError("This structure is on the ring, "
                                    + "so it can't change two substitutions.")
            a_result = ''
        
        frags = addfrag(frags, i, verbose=verbose)
        if substitutions[key] == '':
            pass
        elif side == 1:
            frags[i+1] = '('+ substitutions[key] +')'
            for i0 in range(j):
                del frags[0]
            frags = addfrag(frags, -1)
            frags[0] = substitutions[key]
        elif side == 2:
            frags = addfrag(frags, i, verbose=verbose)
            frags[i+1] = '('+ substitutions[key] +')'
            if keyname == 'ring3' or keyname == 'cyclo5':
                substitutions[key] = 'C{0}CC{0}'.format(big_num+2)
            if keyname == 'ring5' or keyname == 'cyclo5':
                substitutions[key] = 'C{0}CCCC{0}'.format(big_num+2)
            if keyname == 'ring6' or keyname == 'cyclo5':
                substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+2)
            if keyname == 'Benzene' or keyname == 'benzene' or keyname == 'Ph':
                substitutions[key] = 'c{0}ccccc{0}'.format(big_num+2)
            if keyname == 'Naphthalene' or keyname == 'naphthalene' or keyname == 'NA':
                substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+3, big_num+4)
            frags[j+2] = substitutions[key]
    SMILES = merge(frags)
    if key in ('Z', 'E'):
        SMILES = doublebondstyle(SMILES, style, j=i, output='smiles')
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)

def CringSub(substitutions, frags, i, C1, verbose=False, special=False):
    key = list(substitutions.keys())[0]
    if key == 'Z' or key == 'E':
        raise TypeError("This carbon atom is on the C-C bond.")
    elif key == '2':
        raise ValueError("The structure can't be added on this moiety.")
    elif ('R' in keys or 'S' in keys) and i == 0:
        raise TypeError("For the first atom, it can only use '1' or '2'.")
    j = i
    branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j)
    if '@@' in frags[i]:
        property_branch = ['S', 'R']
    elif '@' in frags[i]:
        property_branch = ['R', 'S']
    else:
        property_branch = ['1', '1']
    
    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if verbose == True:
                    print(a.group(0), big_num)
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))

    num = re.search(r'[0-9]+', C1)

    keyname = substitutions[key]
    if substitutions[key] == 'H':
        substitutions[key] = ''
    if substitutions[key] == 'NO2':
        substitutions[key] = '[N+]([O-])=O'
    if substitutions[key] == 'CHO':
        substitutions[key] = 'C(=O)'
    if substitutions[key] == 'COOH':
        substitutions[key] = 'C(=O)O'
    if substitutions[key] == 'CN':
        substitutions[key] = 'C#N'
    if substitutions[key] == 'SO3H':
        substitutions[key] = 'S(=O)(=O)O'
    if substitutions[key] == 'PO3H2':
        substitutions[key] = 'P(=O)(O)O'
    if substitutions[key] == 'CF3':
        substitutions[key] = 'C(F)(F)F'
    if substitutions[key] == 'CCl3':
        substitutions[key] = 'C(Cl)(Cl)Cl'
    if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
        substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
    if substitutions[key] == 'ring5' or substitutions[key] == 'cyclo5':
        substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
    if substitutions[key] == 'ring6' or substitutions[key] == 'cyclo6':
        substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
    if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
        substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
    if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
        substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)

    a2_result = ''
    a = re.match(r'\(+', frags[i])
    a2 = re.search(r'\)+', frags[i])
    if a2:
        a2_result = a2.group(0)
    if a2 and key in ('R', 'S'):
        substitutions['1'] = substitutions[key]
        del substitutions[key]
        key = '1'
    if key == '1':
        if a:
            frags[i] = a.group(0) + C1
            if branch[1] != '':
                for j in range(branch_frag2):
                    del frags[i+branch_frag1+1]
            a_result = a.group(0)
        else:
            frags[i] = C1
            if branch[1] != '':
                for j in range(branch_frag2):
                    del frags[i+branch_frag1+1]
            a_result = ''
        frags = addfrag(frags, i+branch_frag1, verbose=verbose)
        if substitutions[key] == '':
            pass
        elif i+1 == len(frags)-1:
            frags[i+branch_frag1+1] = substitutions[key]
        else:
            frags[i+branch_frag1+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
    elif key == 'R':
        if substitutions[key] == '':
            if a:
                if branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = a.group(0) + '[C@H]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    frags[i] = '[C@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = '[C@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = '[C@H]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = ''

        elif substitutions[key] != '':
            if a:
                if branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = a.group(0) + '[C@]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    frags[i] = '[C@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = '[C@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = '[C@]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            elif i+1 == len(frags)-1:
                frags[i+1] = substitutions[key] + len(a2_result) * ')'
            else:
                frags[i+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
    elif key == 'S':
        if substitutions[key] == '':
            if a:
                if branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = a.group(0) + '[C@@H]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    frags[i] = '[C@@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = '[C@@H]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = '[C@@H]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = ''

        elif substitutions[key] != '':
            if a:
                if branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = a.group(0) + '[C@@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = a.group(0) + '[C@@]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = a.group(0)
            else:
                if branch_frag1 == 0:
                    frags[i] = '[C@@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 == 0:
                    frags[i] = '[C@@]' + num.group(0)
                elif branch[0] != '' and branch_frag1 != 0:
                    frags[i] = '[C@@]' + num.group(0)
                    for j in range(branch_frag1):
                        del frags[i+1]
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            elif i+1 == len(frags)-1:
                frags[i+1] = substitutions[key] + len(a2_result) * ')'
            else:
                frags[i+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
    SMILES = merge(frags)
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)

def NSub(substitutions, frags, i, verbose=False, special=False):
    key = list(substitutions.keys())[0]
    if key == 'Z' or key == 'E':
        raise TypeError("This is the nitrogen atom.")
    elif key == 'R' or key == 'S':
        raise TypeError("This is the nitrogen atom.")
    elif key == '2':
        if '=' in frags[i]:
            raise TypeError("This structure can only add one substitution.")
        try:
            if '=' in frags[i+1]:
                raise TypeError("This structure can only add one substitution.")
        except:
            pass
    
    j = i
    branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j, verbose=verbose)

    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if verbose == True:
                    print(a.group(0), big_num)
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))

    keyname = substitutions[key]
    if substitutions[key] == 'H':
        substitutions[key] = ''
    if substitutions[key] == 'NO2':
        substitutions[key] = '[N+]([O-])=O'
    if substitutions[key] == 'CHO':
        substitutions[key] = 'C(=O)'
    if substitutions[key] == 'COOH':
        substitutions[key] = 'C(=O)O'
    if substitutions[key] == 'CN':
        substitutions[key] = 'C#N'
    if substitutions[key] == 'SO3H':
        substitutions[key] = 'S(=O)(=O)O'
    if substitutions[key] == 'PO3H2':
        substitutions[key] = 'P(=O)(O)O'
    if substitutions[key] == 'CF3':
        substitutions[key] = 'C(F)(F)F'
    if substitutions[key] == 'CCl3':
        substitutions[key] = 'C(Cl)(Cl)Cl'
    if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
        substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
    if substitutions[key] == 'ring5' or substitutions[key] == 'cyclo5':
        substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
    if substitutions[key] == 'ring6' or substitutions[key] == 'cyclo6':
        substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
    if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
        substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
    if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
        substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)

    a2_result = ''
    a = re.match(r'\(+?', frags[i])
    a2 = re.search(r'\)+', frags[i])
    if a2:
        a2_result = a2.group(0)
    if key == '1':
        if frags[i+1][0] == '=' and i != len(frags)-1:
            if a:
                frags[i] = a.group(0) + '[N+]'
                a_result = a.group(0)
            else:
                frags[i] = '[N+]'
                a_result = ''
        elif frags[i][0] == '=':
            if a:
                frags[i] = a.group(0) + '=[N+]'
                a_result = a.group(0)
            else:
                frags[i] = '=[N+]'
                a_result = ''
        elif a:
            frags[i] = a.group(0) + 'N'
            if branch[0] != '':
                for j in range(branch_frag1):
                    del frags[i+1]
            a_result = a.group(0)
        else:
            frags[i] = 'N'
            if branch[0] != '':
                for j in range(branch_frag1):
                    del frags[i+1]
            a_result = ''
        frags = addfrag(frags, i+branch_frag1, verbose=verbose)
        if substitutions[key] == '':
            pass
        elif i+1 == len(frags)-1:
            frags[i+branch_frag1+1] = substitutions[key]
            i = i + 1
        else:
            frags[i+branch_frag1+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
            i = i + 1
    elif key == '2':
        if frags[i][0] == '=':
            if a:
                frags[i] = a.group(0) + '=[N+]'
                a_result = a.group(0)
            else:
                frags[i] = '=[N+]'
                a_result = ''
        elif a:
            if len(a.group(0)) <= len(a2.group(0)):
                frags[i] = a.group(0) + 'N'
                a_result = a.group(0)
            else:
                frags[i] = a.group(0) + '[N+]'
                a_result = a.group(0)
        else:
            if len(a2_result) > 0:
                frags[i] = 'N'
                a_result = ''
            elif len(a2_result) == 0 and i == len(frags) - 1:
                frags[i] = 'N'
                a_result = ''
            else:
                frags[i] = '[N+]'
                a_result = ''
        frags = addfrag(frags, i, verbose=verbose)
        frags = addfrag(frags, i, verbose=verbose)
        if substitutions[key] == '':
            pass
        else:
            frags[i+1] = '('+ substitutions[key] +')'
            if keyname == 'ring3' or keyname == 'cyclo3':
                substitutions[key] = 'C{0}CC{0}'.format(big_num+2)
            if keyname == 'ring5' or keyname == 'cyclo5':
                substitutions[key] = 'C{0}CCCC{0}'.format(big_num+2)
            if keyname == 'ring6' or keyname == 'cyclo6':
                substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+2)
            if keyname == 'Benzene' or keyname == 'benzene' or keyname == 'Ph':
                substitutions[key] = 'c{0}ccccc{0}'.format(big_num+2)
            if keyname == 'Naphthalene' or keyname == 'naphthalene' or keyname == 'NA':
                substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+3, big_num+4)
            if i+2 == len(frags)-1:
                frags[i+2] = substitutions[key] + len(a2_result) * ')'
            else:
                frags[i+2] = '('+ substitutions[key] +')' + len(a2_result) * ')'
    SMILES = merge(frags)
    if verbose == True:
        print(SMILES)
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)

def XtraNSub(substitutions, frags, i, verbose=False, special=False):
    keys = list(substitutions.keys())
    big_num = 0
    if len(keys) == 2 and ({*keys} != {'R', 'S'}, {*keys} != {'1', '1'}):
        raise TypeError("For two substitutions, it can only use 'R' and 'S' or two '1's.")
    elif '1' in keys and '=' not in frags[i]:
        raise TypeError("For this structure, '1' can't show the substitution location.")
    elif '1' in keys and '=' not in frags[i+1] and i != len(frags) - 1:
        raise TypeError("For this structure, '1' can't show the substitution location.")
    
    j = i
    branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j, verbose=verbose)
    if verbose == True:
        print(branch)

    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))
    
    for key in keys:
        if key == 'Z' or key == 'E':
            raise TypeError("This carbon atom is on the C-C bond.")
        keyname = substitutions[key]
        if substitutions[key] == 'H':
            substitutions[key] = ''
        if substitutions[key] == 'NO2':
            substitutions[key] = '[N+]([O-])=O'
        if substitutions[key] == 'CHO':
            substitutions[key] = 'C(=O)'
        if substitutions[key] == 'COOH':
            substitutions[key] = 'C(=O)O'
        if substitutions[key] == 'CN':
            substitutions[key] = 'C#N'
        if substitutions[key] == 'SO3H':
            substitutions[key] = 'S(=O)(=O)O'
        if substitutions[key] == 'PO3H2':
            substitutions[key] = 'P(=O)(O)O'
        if substitutions[key] == 'CF3':
            substitutions[key] = 'C(F)(F)F'
        if substitutions[key] == 'CCl3':
            substitutions[key] = 'C(Cl)(Cl)Cl'
        if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
            substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
        if substitutions[key] == 'ring5' or substitutions[key] == 'cyclo5':
            substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
        if substitutions[key] == 'ring6' or substitutions[key] == 'cyclo6':
            substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
        if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
            substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
        if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
            substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)
        a2_result = ''
        a = re.match(r'\(+', frags[i])
        a2 = re.search(r'\)+', frags[i])
        if a2:
            a2_result = a2.group(0)
        if a2 and key in ('R', 'S'):
            substitutions['1'] = substitutions[key]
            del substitutions[key]
            key = '1'
        if key == '1':
            if a:
                frags[i] = a.group(0) + '[N+]'
                for j in range(branch_frag1):
                    del frags[i+1]
                a_result = a.group(0)
            else:
                frags[i] = '[N+]'
                for j in range(branch_frag1):
                    del frags[i+1]
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            else:
                frags[i+1] = '(' + substitutions[key] +')' + len(a2_result) * ')'
            i = i + 1
        elif key == 'R':
            if a:
                if '' not in branch:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    frags[i] = a.group(0) + '[N+]'
                a_result = a.group(0)
            else:
                if verbose == True:
                    print(branch_frag1)
                if '' not in branch:
                    for j in range(branch_frag1):
                        del frags[i+1]
                    frags[i] = '[N+]'
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            else:
                frags[i+1] = '(' + substitutions[key] + ')' + len(a2_result) * ')'
            i = i + 1
        elif key == 'S':
            if a:
                if '' not in branch:
                    for j in range(branch_frag2):
                        del frags[i+branch_frag1+1]
                    frags[i] = a.group(0) + '[N+]'
                a_result = a.group(0)
            else:
                if verbose == True:
                    print(branch_frag1)
                if '' not in branch:
                    for j in range(branch_frag1):
                        del frags[i+branch_frag1+1]
                    frags[i] = '[N+]'
                a_result = ''
            frags = addfrag(frags, i+branch_frag1, verbose=verbose)
            if substitutions[key] == '':
                pass
            else:
                frags[i+branch_frag1+1] = '(' + substitutions[key] + ')' + len(a2_result) * ')'
            j = j + 1
        elif key == '2':
            if a:
                if '' not in branch:
                    for j in range(branch_frag1+branch_frag2):
                        del frags[i+1]
                    frags[i] = a.group(0) + '[N+]'
                a_result = a.group(0)
            else:
                if '' not in branch:
                    for j in range(branch_frag1+branch_frag2):
                        del frags[i+1]
                    frags[i] = '[N+]'
                a_result = ''
            frags = addfrag(frags, i, verbose=verbose)
            frags = addfrag(frags, i, verbose=verbose)
            if substitutions[key] == '':
                pass
            else:
                frags[i+1] = '('+ substitutions[key] +')'
                if keyname == 'ring3' or keyname == 'cyclo3':
                    substitutions[key] = 'C{0}CC{0}'.format(big_num+2)
                if keyname == 'ring5' or keyname == 'cyclo5':
                    substitutions[key] = 'C{0}CCCC{0}'.format(big_num+2)
                if keyname == 'ring6'or keyname == 'cyclo6':
                    substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+2)
                if keyname == 'Benzene' or keyname == 'benzene' or keyname == 'Ph':
                    substitutions[key] = 'c{0}ccccc{0}'.format(big_num+2)
                if keyname == 'Naphthalene' or keyname == 'naphthalene' or keyname == 'NA':
                    substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+3, big_num+4)
                frags[i+2] = '('+ substitutions[key] +')' + len(a2_result) * ')'
            
    SMILES = merge(frags)
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)

def OSub(substitutions, frags, i, verbose=False, special=False):
    key = list(substitutions.keys())[0]
    if key == 'Z' or key == 'E':
        raise TypeError("This is the oxygen atom.")
    elif key == 'R' or key == 'S':
        raise TypeError("This is the oxygen atom.")
    elif key == '2':
        raise TypeError("This oxygen atom can't add two substitutions.")
    j = i
    
    branch, brack_left, brack_right, branch_frag1, branch_frag2 = bracket_module(frags, j, verbose=verbose)

    big_num = 0
    for index in frags:
        a = re.search(r'\+?\-?[0-9]+\+?\-?', index)
        if a:
            if '+' in a.group(0):
                break
            elif '-' in a.group(0):
                break
            else:
                if verbose == True:
                    print(a.group(0), big_num)
                if eval(a.group(0)) > big_num:
                    big_num = eval(a.group(0))

    keyname = substitutions[key]
    if substitutions[key] == 'H':
        substitutions[key] = ''
    if substitutions[key] == 'NO2':
        substitutions[key] = '[N+]([O-])=O'
    if substitutions[key] == 'CHO':
        substitutions[key] = 'C(=O)'
    if substitutions[key] == 'COOH':
        substitutions[key] = 'C(=O)O'
    if substitutions[key] == 'CN':
        substitutions[key] = 'C#N'
    if substitutions[key] == 'SO3H':
        substitutions[key] = 'S(=O)(=O)O'
    if substitutions[key] == 'PO3H2':
        substitutions[key] = 'P(=O)(O)O'
    if substitutions[key] == 'CF3':
        substitutions[key] = 'C(F)(F)F'
    if substitutions[key] == 'CCl3':
        substitutions[key] = 'C(Cl)(Cl)Cl'
    if substitutions[key] == 'ring3' or substitutions[key] == 'cyclo3':
        substitutions[key] = 'C{0}CC{0}'.format(big_num+1)
    if substitutions[key] == 'ring5' or substitutions[key] == 'cyclo5':
        substitutions[key] = 'C{0}CCCC{0}'.format(big_num+1)
    if substitutions[key] == 'ring6' or substitutions[key] == 'cyclo6':
        substitutions[key] = 'C{0}CCCCC{0}'.format(big_num+1)
    if substitutions[key] == 'Benzene' or substitutions[key] == 'benzene' or substitutions[key] == 'Ph':
        substitutions[key] = 'c{0}ccccc{0}'.format(big_num+1)
    if substitutions[key] == 'Naphthalene' or substitutions[key] == 'naphthalene' or substitutions[key] == 'NA':
        substitutions[key] = 'c{0}cc{1}ccccc{1}cc{0}'.format(big_num+1, big_num+2)

    a2_result = ''
    a = re.match(r'\(+?', frags[i])
    a2 = re.search(r'\)+', frags[i])
    if a2:
        a2_result = a2.group(0)
    if key == '1':
        if a:
            if len(a.group(0)) > len(a2_result):
                raise TypeError("This oxygen atom can't add any substitutions.")
            if 'O' in frags[i]:
                frags[i] = a.group(0) + 'O'
            elif 'S' in frags[i]:
                frags[i] = a.group(0) + 'S'
        else:
            if a2:
                if 'O' in frags[i]:
                    frags[i] = 'O'
                elif 'S' in frags[i]:
                    frags[i] = 'S'
            if not a2 and i != len(frags)-1:
                raise TypeError("This oxygen atom can't add any substitutions.")
        frags = addfrag(frags, i, verbose=verbose)
        if substitutions[key] == '':
            pass
        elif i+1 == len(frags)-1:
            frags[i+branch_frag1+1] = substitutions[key]
            i = i + 1
        else:
            frags[i+branch_frag1+1] = '(' +substitutions[key] +')' + len(a2_result) * ')'
            i = i + 1
        #else:
        #    frags[i+branch_frag1+1] =  substitutions[key] + len(a2_result) * ')'
        #    i = i + 1
    SMILES = merge(frags)
    if special == True:
        return SMILES
    else:
        mole_SMILES = Chem.MolFromSmiles(SMILES)
        #display(mole_SMILES)
        return Chem.MolToSmiles(mole_SMILES)