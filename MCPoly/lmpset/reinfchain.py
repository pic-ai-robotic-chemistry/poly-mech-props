import re
import os
import sys
import warnings
import math as m
import numpy as np
from ase.io import read,  write

mydir = os.path.dirname( __file__ )
lmpdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(lmpdir)
from infchain import infchain
from DataToXyz import DataToXyz

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'status')
sys.path.append(statusdir)
from molesplit import molesplit

def ch1(num):
    if num >= 0:
        return num + 1
    else:
        return num - 1

def near(atomnum, atoms, tolerance=2.000):
    atoms_not_H = []
    for i, atom in enumerate(atoms):
        if atom.symbol != 'H':
            atoms_not_H.append(i)
    mole = [atomnum]
    pre_distances = atoms.get_distances(atomnum, atoms_not_H)
    distances = np.delete(pre_distances, np.where(pre_distances >= tolerance))
    distances = np.delete(distances, np.where(distances == 0))
    distances.sort()
    length = len(distances)
    if length > 4:
        distances = distances[:4]
    for i, atom in enumerate(atoms_not_H):
        distance = atoms.get_distance(atom, atomnum)
        if distance in distances:
            mole.append(atom)
    return mole

def near_between(atomnum, atoms, maxi=2.000, mini=1.330):
    atoms_ = []
    for i, atom in enumerate(atoms):
        atoms_.append(i)
    mole = [atomnum]
    pre_distances = atoms.get_distances(atomnum, atoms_)
    distances = np.delete(pre_distances, np.where(pre_distances >= maxi))
    distances = np.delete(distances, np.where(distances <= mini))
    distances.sort()
    length = len(distances)
    if length > 4:
        distances = distances[:4]
    for i, atom in enumerate(atoms_):
        distance = atoms.get_distance(atom, atomnum)
        if distance in distances:
            mole.append(atom)
    return mole

def find_midBADI(file, loc='./'):
    opath = os.getcwd()
    os.chdir(loc)
    DataToXyz(file)
    atoms = read(file+'.xyz')
    for i,atom in enumerate(atoms):
        if atom.symbol != 'H':
            break
    mole = molesplit(file, atomnum=i)
    atom = mole[-1]
    mole = molesplit(file, atomnum=atom)
    pos_a = atoms[mole[0]].position
    pos_b = atoms[mole[1]].position
    if pos_a[0] < pos_b[0]:
        mole.reverse()
    atomnum1 = mole[len(mole)//2 - 1]
    atomnum2 = mole[len(mole)//2]
    index1 = near(atomnum1, atoms)
    index2 = near(atomnum2, atoms)
    i = 2
    head = mole[0]
    tail = mole[-1]
    while 1:
        if pos_a[0] < pos_b[0]:
            if atomnum1 not in index2:
                atomnum2 = mole[len(mole)//2 + i - 1]
                i = i + 1
                index2 = near(atomnum2, atoms)
            else:
                break
        else:
            if atomnum2 not in index1:
                atomnum1 = mole[len(mole)//2 - i]
                i = i + 1
                index1 = near(atomnum1, atoms)
            else:
                break
    i = 1
    index1 = near(head, atoms)
    index2 = near(tail, atoms)
    index1_2 = near_between(head, atoms, mini=0)
    index2_2 = near_between(tail, atoms, mini=0)
    while 1:
        if pos_a[0] < pos_b[0]:
            if len(index1) == 2 and len(index1_2) < 5:
                head = mole[i]
                i = i + 1
                index1 = near(head, atoms)
                index1_2 = near_between(head, atoms, mini=0)
            else:
                break
        else:
            if len(index2) == 2 and len(index2_2) < 5:
                tail = mole[-i-1]
                i = i + 1
                index2 = near(tail, atoms)
                index2_2 = near_between(tail, atoms, mini=0)
            else:
                break
    pos_z = atoms[-1].position
    pos_c = atoms[atomnum1].position
    pos_d = atoms[atomnum2].position
    atomnum1 = atomnum1 + 1
    atomnum2 = atomnum2 + 1
    #if pos_a[0] < pos_b[0]:
    xxx = head
    head = tail
    tail = xxx
    
    print('MIDDLE', atomnum1, atomnum2)
    print('RESULT', head, tail)
    os.system('rm {0}.xyz'.format(file))
    
    f = open('{0}.data'.format(file), 'r')
    b1 = 0
    b2 = 0
    b3 = 0
    b4 = 0
    bonds_a = []
    angles_a = []
    dihedrals_a = []
    impropers_a = []
    bonds_b = []
    angles_b = []
    dihedrals_b = []
    impropers_b = []
    
    for line in f:
        a1 = re.search('Bonds', line)
        a2 = re.search('Angles', line)
        a3 = re.search('Dihedrals', line)
        a4 = re.search('Impropers', line)
        if a1:
            b1 = b1 + 1
        elif b1 == 1:
            b1 = b1 + 1
        elif b1 == 2:
            index = re.findall(r'[0-9]+', line)
            if str(atomnum1) in index[2:]:
                if str(atomnum2) in index[2:]:
                    bonds_a.append(index[1])
                    bonds_b.append(index[1])
            if index == []:
                b1 = 0
                continue
        if a2:
            b2 = b2 + 1
        elif b2 == 1:
            b2 = b2 + 1
        elif b2 == 2:
            index = re.findall(r'[0-9]+', line)
            if str(atomnum1) in index[2:]:
                if str(atomnum2) in index[2:]:
                    angles_a.append([index[1], atoms[eval(index[2])-1].symbol,\
                                     atoms[eval(index[3])-1].symbol, atoms[eval(index[4])-1].symbol])
                    angles_b.append(index[1:])
            if index == []:
                b2 = 0
                continue
        if a3:
            b3 = b3 + 1
        elif b3 == 1:
            b3 = b3 + 1
        elif b3 == 2:
            index = re.findall(r'[0-9]+', line)
            if str(atomnum1) in index[2:]:
                if str(atomnum2) in index[2:]:
                    dihedrals_a.append([index[1], atoms[eval(index[2])-1].symbol,\
                                        atoms[eval(index[3])-1].symbol, atoms[eval(index[4])-1].symbol,\
                                        atoms[eval(index[5])-1].symbol])
                    dihedrals_b.append(index[1:])
            if index == []:
                b3 = 0
                continue
        if a4:
            b4 = b4 + 1
        elif b4 == 1:
            b4 = b4 + 1
        elif b4 == 2:
            index = re.findall(r'[0-9]+', line)
            if str(atomnum1) in index[2:]:
                if str(atomnum2) in index[2:]:
                    impropers_a.append([index[1], atoms[eval(index[2])-1].symbol,\
                                        atoms[eval(index[3])-1].symbol, atoms[eval(index[4])-1].symbol,\
                                        atoms[eval(index[5])-1].symbol])
                    impropers_b.append(index[1:])
            if index == []:
                b4 = 0
                continue
    f.close()
        
    for i,angle in enumerate(angles_a):
        for j,ele in enumerate(angle):
            if j == 0:
                continue
            if angles_b[i][j] == str(atomnum1):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    angles_b[i][j] = str(atomnum2)
                else:
                    pass
            elif angles_b[i][j] == str(atomnum2):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    angles_b[i][j] = str(atomnum1)
                else:
                    pass
            else:
                angles_b[i][j] = angles_a[i][j]
    for i,dihedral in enumerate(dihedrals_a):
        for j,ele in enumerate(dihedral):
            if j == 0:
                continue
            if dihedrals_b[i][j] == str(atomnum1):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    dihedrals_b[i][j] = str(atomnum2)
                else:
                    pass
            elif dihedrals_b[i][j] == str(atomnum2):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    dihedrals_b[i][j] = str(atomnum1)
                else:
                    pass
            else:
                dihedrals_b[i][j] = dihedrals_a[i][j]
    for i,improper in enumerate(impropers_a):
        for j,ele in enumerate(improper):
            if j == 0:
                continue
            if impropers_b[i][j] == str(atomnum1):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    impropers_b[i][j] = str(atomnum2)
                else:
                    pass
            elif impropers_b[i][j] == str(atomnum2):
                if (atomnum1-atomnum2) * (head-tail) < 0:
                    impropers_b[i][j] = str(atomnum1)
                else:
                    pass
            else:
                impropers_b[i][j] = impropers_a[i][j]
    os.chdir(opath)
    return [head+1, tail+1, bonds_b, angles_b, dihedrals_b, impropers_b]

def find_angle(file, head, tail, angles, bonds, xminnum,\
               xmaxnum, realxminnum, realxmaxnum):
    DataToXyz(file)
    atoms = read(file+'.xyz')
    results = []
    if head - tail > 0:
        x = head
        head = tail
        tail = x
    if realxminnum - realxmaxnum > 0:
        x = realxminnum
        realxminnum = realxmaxnum
        realxmaxnum = x
    for bond in bonds:
        for angle in angles:
            if {head, tail} == {eval(bond[2]), eval(bond[3])}:
                continue
            elif str(xminnum) in bond[2:]:
                continue
            elif str(xmaxnum) in bond[2:]:
                continue
            elif head == eval(bond[2]):
                #print(head, tail, bond[2:], angle, 1)
                t1 = eval(bond[3])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) < eval(angle[3]):
                    results.append([angle[0], t1, realxminnum, realxmaxnum])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) < eval(angle[1]):
                    results.append([angle[0], t1, realxminnum, realxmaxnum])
                    break
            elif head == eval(bond[3]):
                #print(head, tail, bond[2:], angle, 2)
                t1 = eval(bond[2])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) < eval(angle[3]):
                    results.append([angle[0], realxminnum, t1, realxmaxnum])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) < eval(angle[1]):
                    results.append([angle[0], realxmaxnum, realxminnum, t1])
                    break
            elif tail == eval(bond[2]):
                #print(head, tail, bond[2:], angle, 3)
                t1 = eval(bond[3])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) > eval(angle[3]):
                    results.append([angle[0], t1, realxmaxnum, realxminnum])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) > eval(angle[1]):
                    results.append([angle[0], realxminnum, realxmaxnum, t1])
                    break
            elif tail == eval(bond[3]):
                #print(head, tail, bond[2:], angle, 4)
                t1 = eval(bond[2])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) > eval(angle[3]):
                    results.append([angle[0], t, realxminnum, realxmaxnum])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) > eval(angle[1]):
                    results.append([angle[0], realxminnum, realxmaxnum, t1])
                    break
    os.system('rm {0}.xyz'.format(file))
    print(results)
    return results

def find_dihedral(file, head, tail, dihedrals, extra_angles,\
                  angles, xminnum, xmaxnum, realxminnum, realxmaxnum):
    DataToXyz(file)
    atoms = read(file+'.xyz')
    results = []
    if head - tail > 0:
        x = head
        head = tail
        tail = x
    if realxminnum - realxmaxnum > 0:
        x = realxminnum
        realxminnum = realxmaxnum
        realxmaxnum = x
    for angle in extra_angles:
        for angle2 in extra_angles:
            if angle == angle2:
                continue
            for dihedral in dihedrals:
                if angle[1:3] == angle2[2:]:
                    t1 = ch1(angle2[1])
                    t2 = ch1(angle[3])
                    if atoms[t1-1].symbol == dihedral[1] and atoms[t2-1].symbol == dihedral[4]:
                        #print(realxminnum - num2, realxmaxnum - num, angle2[2], angle2[3])
                        if 1:#(realxminnum - num2, realxmaxnum - num) == (angle2[2], angle2[3]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle2[1], angle2[2], angle2[3], angle[3]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                            break
                    elif atoms[t1-1].symbol == dihedral[4] and atoms[t2-1].symbol == dihedral[1]:
                        #print(realxminnum - num2, realxmaxnum - num, angle2[2], angle2[3])
                        if 1:# (realxminnum - num2, realxmaxnum - num) == (angle2[3], angle2[2]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle2[1], angle2[2], angle2[3], angle[3]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                            break
                elif angle[1:3] == angle2[2:0:-1]:
                    t1 = ch1(angle2[3])
                    t2 = ch1(angle[3])
                    if atoms[t1-1].symbol == dihedral[1] and atoms[t2-1].symbol == dihedral[4]:
                        #print(realxminnum - num2, realxmaxnum - num, angle[1], angle2[2])
                        if 1:# (realxminnum - num2, realxmaxnum - num) == (angle[1], angle[2]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle2[3], angle[1], angle[2], angle[3]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                        break
                    elif atoms[t1-1].symbol == dihedral[4] and atoms[t2-1].symbol == dihedral[1]:
                        #print(realxminnum - num2, realxmaxnum - num, angle[1], angle[2])
                        if 1:# (realxminnum - num2, realxmaxnum - num) == (angle[2], angle[1]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle2[3], angle[1], angle[2], angle[3]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                        break
                elif angle[2:] == angle2[3:1:-1]:
                    t1 = ch1(angle[1])
                    t2 = ch1(angle2[1])
                    if atoms[t1-1].symbol == dihedral[1] and atoms[t2-1].symbol == dihedral[4]:
                        #print(realxminnum - num2, realxmaxnum - num, angle[2], angle[3])
                        if 1:# (realxminnum - num2, realxmaxnum - num) == (angle[2], angle[3]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle[1], angle[2], angle[3], angle2[1]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                        break
                    elif atoms[t1-1].symbol == dihedral[4] and atoms[t2-1].symbol == dihedral[1]:
                        #print(realxminnum - num2, realxmaxnum - num, angle[2], angle[3])
                        if 1:# (realxminnum - num2, realxmaxnum - num) == (angle[3], angle[2]):
                            #print(t1, t2, dihedral, angle, angle2)
                            result = [dihedral[0], angle[1], angle[2], angle[3], angle2[1]]
                            if result in results or [result[0], *result[-1:0:-1]] in results:
                                pass
                            else:
                                results.append(result)
                        break
    for angle in angles:
        if head in angle and tail in angle:
            continue
        elif str(xminnum) in angle[2:]:
            continue
        elif str(xmaxnum) in angle[2:]:
            continue
        for dihedral in dihedrals:
            if head == eval(angle[2]):
                #print(head, tail, angle[2:], dihedral,1)
                t1 = eval(angle[3])
                t2 = eval(angle[4])
                if atoms[t2-1].symbol == dihedral[1] and atoms[t1-1].symbol == dihedral[2]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t2, t1,\
                                    realxminnum, realxmaxnum])
                    break
                elif atoms[t1-1].symbol == dihedral[3] and atoms[t2-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t2, t1,\
                                    realxminnum, realxmaxnum])
                    break
            if head == eval(angle[4]):
                #print(head, tail, angle[2:], dihedral,2)
                t1 = eval(angle[2])
                t2 = eval(angle[3])
                if atoms[t1-1].symbol == dihedral[1] and atoms[t2-1].symbol == dihedral[2]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t1, t2,\
                                    realxminnum, realxmaxnum])
                    break
                elif atoms[t1-1].symbol == dihedral[4] and atoms[t2-1].symbol == dihedral[3]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t1, t2,\
                                    realxminnum, realxmaxnum])
                    break
            if tail == eval(angle[2]):
                #print(head, tail, angle[2:], dihedral,3)
                t1 = eval(angle[3])
                t2 = eval(angle[4])
                if atoms[t1-1].symbol == dihedral[3] and atoms[t2-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum, realxmaxnum,\
                                    t1, t2])
                    break
                elif atoms[t1-1].symbol == dihedral[2] and atoms[t2-1].symbol == dihedral[1]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum, realxmaxnum,\
                                    t1, t2])
                    break
            if tail == eval(angle[4]):
                #print(head, tail, angle[2:], dihedral,4)
                t1 = eval(angle[2])
                t2 = eval(angle[3])
                if atoms[t2-1].symbol == dihedral[3] and atoms[t1-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum, realxmaxnum,\
                                    t2, t1])
                    break
                elif atoms[t2-1].symbol == dihedral[2] and atoms[t1-1].symbol == dihedral[1]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum, realxmaxnum,\
                                    t2, t1])
                    break
    os.system('rm {0}.xyz'.format(file))
    print(results)
    return results

def find_improper(file, head, tail, impropers, extra_angles,\
                  angles, xminnum, xmaxnum, realxminnum, realxmaxnum):
    DataToXyz(file)
    atoms = read(file+'.xyz')
    num = len(atoms)
    results = []
    if head - tail > 0:
        x = head
        head = tail
        tail = x
    if realxminnum - realxmaxnum > 0:
        x = realxminnum
        realxminnum = realxmaxnum
        realxmaxnum = x
    has = 0
    if atoms[head-1].symbol == 'C':
        haslimit = 2
    if atoms[head-1].symbol == 'N':
        haslimit = 1
    if atoms[head-1].symbol in ('O', 'S'):
        haslimit = 0
    ambh = head
    ambt = tail
    if head in extra_angles[0]:
        head = head
    elif head - num in extra_angles[0]:
        head = head - num
    if tail in extra_angles[0]:
        tail = tail
    elif tail - num in extra_angles[0]:
        tail = tail - num
    for angle in extra_angles:
        for angle2 in extra_angles:
            if angle == angle2:
                continue
            for improper in impropers:
                if has == haslimit:
                    break
                if angle[2] == angle2[2] == head:
                    t1 = angle[1]
                    t2 = angle[3]
                    t3 = angle2[1]
                    t4 = angle2[3]
                    t = [t1, t2, t3, t4]
                    print(angle[1:], angle2[1:], t, improper, [head,tail])
                    t.remove(tail)
                    t.remove(tail)
                    if (atoms[t[0]].symbol in improper[2:] and
                        atoms[t[1]].symbol in improper[2:]):
                        results.append([improper[0], head, tail, t[0], t[1]])
                        #print(results[-1])
                        has = has + 1
                        break
    has = 0
    if atoms[tail-1].symbol == 'C':
        haslimit = 2
    if atoms[tail-1].symbol == 'N':
        haslimit = 1
    if atoms[tail-1].symbol in ('O', 'S'):
        haslimit = 0
        print(haslimit)
    for angle in extra_angles:
        for angle2 in extra_angles:
            if angle == angle2:
                continue
            for improper in impropers:
                if has == haslimit:
                    break
                if angle[2] == angle2[2] == tail:
                    t1 = angle[1]
                    t2 = angle[3]
                    t3 = angle2[1]
                    t4 = angle2[3]
                    t = [t1, t2, t3, t4]
                    print(angle[1:],angle2[1:],t, improper, [head,tail])
                    t.remove(head)
                    t.remove(head)
                    if (atoms[t[0]].symbol in improper[2:] and
                        atoms[t[1]].symbol in improper[2:]):
                        results.append([improper[0], tail, head, t[0], t[1]])
                        has = has + 1
                        #print(results[-1])
                        break
    os.system('rm {0}.xyz'.format(file))
    return results

def reinfchain(file, loc='./', savename=''):
    """
    The method to create a single molecule for periodical chain, specialised for polymers.
    reinfchain(file, bondtype, loc = './', savename='')
    file: Your molecule system name on your .data file.
    bondtype: The bond type between the start and the end of the polymer. 
    loc: File Location. The default is your current location.
    savename: The name of the created file. The default is with '_Chain'
    Example:
        Input:
            from MCPoly.lmpset import chain
            infchain('Poly1', 11)
        Output in Poly1_Chain.data:
            82 atoms
            84 bonds
            154 angles
            211 dihedrals
            47 impropers
            
            9 atom types
            12 bond types
            21 angle types
            31 dihedral types
            3 improper types
            
            ...
            84 11 79 1
            
            ...
    """
    opath = os.getcwd()
    os.chdir(loc)
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    n1 = 0
    t2 = 0
    t3 = 0
    t4 = 0
    t5 = 0
    H = []
    lx = []
    lost2_u = []
    lost2_d = []
    lost3_u = []
    lost3_d = []
    lost4_u = []
    lost4_d = []
    lost5_u = []
    lost5_d = []
    Hloc = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    l = []
    real = [0, 0, 0]
    stat = find_midBADI(file, loc='./')
    print(stat)
    head = stat[0]
    tail = stat[1]
    bonds = stat[2]
    angles = stat[3]
    dihedrals = stat[4]
    impropers = stat[5]
    inp = open('{0}.data'.format(file), 'r')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            num2 = b2.group(0)
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            num5 = b5.group(0)
        if n1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                continue
            else:
                n1 = 2
        if n1 == 2:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if l1 == []:
                n1 = 0
                continue
            if eval(l1[-1]) < 1.1:
                H.append(eval(l1[0]))
        if m1:
            n1 = 1
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if c1:
            w1 = 1
        elif c2:
            w2 = 1
        elif c3:
            w3 = 1
        elif c4:
            w4 = 1
        elif c5:
            w5 = 1
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if l1[0] == num:
                    w1 = 2
        elif w1 == 2:
            xmin = 99999
            xmax = -99999
            xmaxnum = 0
            xminnum = 0
            realxmin = 99999
            realxmax = -99999
            realxmaxnum = 0
            realxminnum = 0
            for i4 in range(len(l)):
                if eval(l[i4][2]) not in H:
                    continue
                if eval(l[i4][-3]) > xmax:
                    xmax = eval(l[i4][-3])
                    xmaxnum = eval(l[i4][0])
                if eval(l[i4][-3]) < xmin:
                    xmin = eval(l[i4][-3])
                    xminnum = eval(l[i4][0])
            for i4 in range(len(l)):
                if eval(l[i4][2]) in H:
                    continue
                if eval(l[i4][-3]) > realxmax:
                    realxmax = eval(l[i4][-3])
                    realxmaxnum = eval(l[i4][0])
                if eval(l[i4][-3]) < realxmin:
                    realxmin = eval(l[i4][-3])
                    realxminnum = eval(l[i4][0])
            w1 = 0
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if l1[0] == num2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t2 = j
            j = 1
            w2 = 0
            k = 0
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if l1[0] == num3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t3 = j
            j = 1
            k = 0
            w3 = 0
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t4 = j
            j = 1
            k = 0
            w4 = 0
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                j = j + 1
            k = k + 1
            l = []
            t5 = j
            j = 1
            k = 0
            w4 = 0
    inp.close()
    #print(t2, t3, t4, t5)
    inp = open('{0}.data'.format(file), 'r')
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    try:
        if savename == '':
            f = open('{0}_Chain.data'.format(file), 'x')
        else:
            f = open('{0}.data'.format(savename), 'x')
    except:
        if savename == '':
            f = open('{0}_Chain.data'.format(file), 'w')
        else:
            f = open('{0}.data'.format(savename), 'w')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        a7 = re.search(r' ylo yhi', line)
        a8 = re.search(r' zlo zhi', line)
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            f.write('{0} atoms\n'.format(eval(b1.group(0))-2))
            num = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            f.write('{0} bonds\n'.format(eval(b2.group(0))-1))
            num2 = b2.group(0)
        elif a3:
            f.write('{0} angles\n'.format(eval(b3.group(0))))
            num3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            f.write('{0} dihedrals\n'.format(eval(b4.group(0))+9))
            num4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            f.write('{0} impropers\n'.format(eval(b5.group(0))))
            num5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1]) - eval(x1[0])
            f.write('{0:>10.5f} {1:>10.5f} xlo xhi\n'.format(xmin, xmax))
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1])-eval(y1[0])
            f.write('{0:>10.5f} {1:>10.5f} ylo yhi\n'.format(eval(y1[0]), eval(y1[1])))
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1]) - eval(z1[0])
            f.write('{0:>10.5f} {1:>10.5f} zlo zhi\n'.format(eval(z1[0]), eval(z1[1])))
        elif c1:
            w1 = 1
            f.write(line)
            f.write('\n')
        elif c2:
            w2 = 1
            f.write(line)
            f.write('\n')
        elif c3:
            w3 = 1
            f.write(line)
            f.write('\n')
        elif c4:
            w4 = 1
            f.write(line)
            f.write('\n')
        elif c5:
            w5 = 1
            f.write(line)
            f.write('\n')
        elif w1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) >= 6:
                l.append(l1)
                if l1[0] == num:
                    w1 = 2
                    j = 1
                    num = eval(num)
        elif w1 == 2:
            for i4 in range(eval(b1.group(0))):
                if j in [xminnum, xmaxnum]:
                    gap = gap + 1
                    j = j + 1
                    continue
                f.write('{0} {1} {2} {3:>7.4f} {4:>10.5f} {5:>10.5f} {6:>10.5f}'.format(j-gap,\
                                                                               eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),\
                                                                               eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])))
                f.write('\n')
                j = j + 1
            l = []
            j = 1
            w1 = 0
            f.write('\n')
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if l1[0] == num2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i4 in range(eval(b2.group(0))):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                   eval(l[i4][2]), eval(l[i4][3])))
                j = j + 1
            k = k + 1
            all_bonds = l
            l = []
            if realxmaxnum in [xminnum, xmaxnum]:
                continue
            elif realxmaxnum > max([xminnum, xmaxnum]):
                realxmaxnum = realxmaxnum - 2
            elif realxmaxnum > min([xminnum, xmaxnum]):
                realxmaxnum = realxmaxnum - 1
            if realxminnum in [xminnum, xmaxnum]:
                continue
            elif realxminnum > max([xminnum, xmaxnum]):
                realxminnum = realxminnum - 2
            elif realxminnum > min([xminnum, xmaxnum]):
                realxminnum = realxminnum - 1
            f.write('{0} {1} {2} {3}\n'.format(j, bonds[0], realxmaxnum, realxminnum))
            j = 1
            w2 = 0
            k = 0
            f.write('\n')
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if l1[0] == num3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                       eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])))
                j = j + 1
            extra_angles = find_angle(file, head, tail, angles, all_bonds,\
                                      xminnum, xmaxnum, realxminnum, realxmaxnum)
            for extra in extra_angles:
                rami_box = [*extra]
                if extra[1] < 0 and extra[1] + num > xmaxnum:
                    extra[1] = extra[1] - 1
                if extra[2] < 0 and extra[2] + num > xmaxnum:
                    extra[2] = extra[2] - 1
                if extra[3] < 0 and extra[3] + num > xmaxnum:
                    extra[3] = extra[3] - 1
                if (extra[1] not in [realxminnum, realxmaxnum]
                        and 0 < extra[1] < min([xminnum, xmaxnum])):
                    extra[1] = extra[1] + 1
                if (extra[2] not in [realxminnum, realxmaxnum]
                        and 0 < extra[2] < min([xminnum, xmaxnum])):
                    extra[2] = extra[2] + 1
                if (extra[3] not in [realxminnum, realxmaxnum]
                        and 0 < extra[3] < min([xminnum, xmaxnum])):
                    extra[3] = extra[3] + 1
                if extra[1] >= max([xminnum, xmaxnum]):
                    extra[1] = extra[1] - 1
                if extra[2] >= max([xminnum, xmaxnum]):
                    extra[2] = extra[2] - 1
                if extra[3] >= max([xminnum, xmaxnum]):
                    extra[3] = extra[3] - 1
                if extra[1] < extra[3]:
                    f.write('{0} {1} {2} {3} {4}\n'.format(j, extra[0], extra[1], extra[2], extra[3]))
                else:
                    f.write('{0} {1} {4} {3} {2}\n'.format(j, extra[0], extra[1], extra[2], extra[3]))
                j = j + 1
            k = k + 1
            all_angles = l
            del all_bonds
            l = []
            j = 1
            k = 0
            w3 = 0
            f.write('\n')
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    continue
                elif eval(l[i4][5]) > max([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-2)
                elif eval(l[i4][5]) > min([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-1)
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                           eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4]), eval(l[i4][5])))
                j = j + 1
            extra_dihedrals = find_dihedral(file, head, tail, dihedrals, extra_angles,\
                                            all_angles, xminnum, xmaxnum, realxminnum, realxmaxnum)
            for extra in extra_dihedrals:
                if extra[1] < 0 and extra[1] + num > xmaxnum:
                    extra[1] = extra[1] - 1
                if extra[2] < 0 and extra[2] + num > xmaxnum:
                    extra[2] = extra[2] - 1
                if extra[3] < 0 and extra[3] + num > xmaxnum:
                    extra[3] = extra[3] - 1
                if extra[4] < 0 and extra[4] + num > xmaxnum:
                    extra[4] = extra[4] - 1
                #if (extra[1] not in [realxminnum, realxmaxnum] and 
                if 0 < extra[1] < min([xminnum, xmaxnum]):
                    extra[1] = extra[1] + 1
                #if (extra[2] not in [realxminnum, realxmaxnum] and
                if 0 < extra[2] < min([xminnum, xmaxnum]):
                    extra[2] = extra[2] + 1
                #if (extra[3] not in [realxminnum, realxmaxnum] and
                if 0 < extra[3] < min([xminnum, xmaxnum]):
                    extra[3] = extra[3] + 1
                #if (extra[4] not in [realxminnum, realxmaxnum] and
                if 0 < extra[4] < min([xminnum, xmaxnum]):
                    extra[4] = extra[4] + 1
                if extra[1] >= max([xminnum, xmaxnum]):
                    extra[1] = extra[1] - 1
                if extra[2] >= max([xminnum, xmaxnum]):
                    extra[2] = extra[2] - 1
                if extra[3] >= max([xminnum, xmaxnum]):
                    extra[3] = extra[3] - 1
                if extra[4] >= max([xminnum, xmaxnum]):
                    extra[4] = extra[4] - 1
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1],\
                                                           extra[2], extra[3], extra[4]))
                j = j + 1
            k = k + 1
            l = []
            all_dihedrals = l
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if l1[0] == num5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            ignh = 0
            ignt = 0
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        ignh = ignh + 1
                    elif eval(l[i4][2]) == xmaxnum:
                        ignt = ignt + 1
                    continue
                elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-2)
                elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                    l[i4][2] = str(eval(l[i4][2])-1)
                if eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        ignh = ignh + 1
                    elif eval(l[i4][3]) == xmaxnum:
                        ignt = ignt + 1
                    continue
                elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-2)
                elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                    l[i4][3] = str(eval(l[i4][3])-1)
                if eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        ignh = ignh + 1
                    elif eval(l[i4][4]) == xmaxnum:
                        ignt = ignt + 1
                    continue
                elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-2)
                elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                    l[i4][4] = str(eval(l[i4][4])-1)
                if eval(l[i4][5]) in [xminnum, xmaxnum]:
                    if eval(l[i4][5]) == xminnum:
                        ignh = ignh + 1
                    elif eval(l[i4][5]) == xmaxnum:
                        ignt = ignt + 1
                    continue
                elif eval(l[i4][5]) > max([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-2)
                elif eval(l[i4][5]) > min([xminnum, xmaxnum]):
                    l[i4][5] = str(eval(l[i4][5])-1)
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                           eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4]), eval(l[i4][5])))
                j = j + 1
            extra_impropers = find_improper(file, head, tail, impropers,\
                                            extra_angles, all_angles,\
                                            xminnum, xmaxnum, realxminnum, realxmaxnum)
            for z,extra in enumerate(extra_impropers):
                if extra[1] < 0 and extra[1] + num > xmaxnum:
                    extra[1] = extra[1] - 1
                if extra[2] < 0 and extra[2] + num > xmaxnum:
                    extra[2] = extra[2] - 1
                if extra[3] < 0 and extra[3] + num > xmaxnum:
                    extra[3] = extra[3] - 1
                if extra[4] < 0 and extra[4] + num > xmaxnum:
                    extra[4] = extra[4] - 1
                #if (extra[1] not in [realxminnum, realxmaxnum] and 
                if 0 < extra[1] < min([xminnum, xmaxnum]):
                    extra[1] = extra[1] + 1
                #if (extra[2] not in [realxminnum, realxmaxnum] and
                if 0 < extra[2] < min([xminnum, xmaxnum]):
                    extra[2] = extra[2] + 1
                #if (extra[3] not in [realxminnum, realxmaxnum] and
                if 0 < extra[3] < min([xminnum, xmaxnum]):
                    extra[3] = extra[3] + 1
                #if (extra[4] not in [realxminnum, realxmaxnum] and
                if 0 < extra[4] < min([xminnum, xmaxnum]):
                    extra[4] = extra[4] + 1
                if extra[1] >= max([xminnum, xmaxnum]):
                    extra[1] = extra[1] - 1
                if extra[2] >= max([xminnum, xmaxnum]):
                    extra[2] = extra[2] - 1
                if extra[3] >= max([xminnum, xmaxnum]):
                    extra[3] = extra[3] - 1
                if extra[4] >= max([xminnum, xmaxnum]):
                    extra[4] = extra[4] - 1
                if z < len(extra_impropers) // 2:
                    if z < len(extra_impropers) // 2 - ignh:
                        continue
                    f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1],\
                                                               extra[2], extra[3], extra[4]))
                else:
                    if z < len(extra_impropers) - ignt:
                        continue
                    f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1],\
                                                               extra[2], extra[3], extra[4]))
                j = j + 1
            k = k + 1
            l = []
            j = 1
            k = 0
            w5 = 0
            f.write('\n')
        else:
            f.write(line)
    f.close()
    inp.close()
    os.chdir(opath)