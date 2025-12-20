import re
import os
import sys
import numpy as np
from ase.io import read

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
    num = len(atoms) - 1
    num2 = len(atoms) - 1
    results = []
    if head - tail > 0:
        x = head
        head = tail
        tail = x
    if realxminnum - realxmaxnum > 0:
        num = 1
    else:
        num2 = 1
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
                    results.append([angle[0], t1 - num2, realxminnum - num2, realxmaxnum - num])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) < eval(angle[1]):
                    results.append([angle[0], t1 - num2, realxminnum - num2, realxmaxnum - num])
                    break
            elif head == eval(bond[3]):
                #print(head, tail, bond[2:], angle, 2)
                t1 = eval(bond[2])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) < eval(angle[3]):
                    results.append([angle[0], realxminnum - num2, t1 - num2, realxmaxnum - num])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) < eval(angle[1]):
                    results.append([angle[0], realxmaxnum - num, realxminnum - num2, t1 - num2])
                    break
            elif tail == eval(bond[2]):
                #print(head, tail, bond[2:], angle, 3)
                t1 = eval(bond[3])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) > eval(angle[3]):
                    results.append([angle[0], t1 - num, realxmaxnum - num, realxminnum - num2])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) > eval(angle[1]):
                    results.append([angle[0], realxminnum - num2, realxmaxnum - num, t1 - num])
                    break
            elif tail == eval(bond[3]):
                #print(head, tail, bond[2:], angle, 4)
                t1 = eval(bond[2])
                if atoms[t1-1].symbol == angle[1] and eval(angle[2]) > eval(angle[3]):
                    results.append([angle[0], t1 - num, realxminnum - num2, realxmaxnum - num])
                    break
                elif atoms[t1-1].symbol == angle[3] and eval(angle[2]) > eval(angle[1]):
                    results.append([angle[0], realxminnum - num2, realxmaxnum - num, t1 - num])
                    break
    os.system('rm {0}.xyz'.format(file))
    print(results)
    return results

def find_dihedral(file, head, tail, dihedrals, extra_angles,\
                  angles, xminnum, xmaxnum, realxminnum, realxmaxnum):
    DataToXyz(file)
    atoms = read(file+'.xyz')
    num = len(atoms) - 1
    num2 = len(atoms) - 1
    results = []
    if head - tail > 0:
        x = head
        head = tail
        tail = x
    if realxminnum - realxmaxnum > 0:
        num = 1
    else:
        num2 = 1
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
                    results.append([dihedral[0], t2 - num2, t1 - num2,\
                                    realxminnum - num2, realxmaxnum - num])
                    break
                elif atoms[t1-1].symbol == dihedral[3] and atoms[t2-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t2 - num2, t1 - num2,\
                                    realxminnum - num2, realxmaxnum - num])
                    break
            if head == eval(angle[4]):
                #print(head, tail, angle[2:], dihedral,2)
                t1 = eval(angle[2])
                t2 = eval(angle[3])
                if atoms[t1-1].symbol == dihedral[1] and atoms[t2-1].symbol == dihedral[2]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t1 - num2, t2 - num2,\
                                    realxminnum - num2, realxmaxnum - num])
                    break
                elif atoms[t1-1].symbol == dihedral[4] and atoms[t2-1].symbol == dihedral[3]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], t1 - num2, t2 - num2,\
                                    realxminnum - num2, realxmaxnum - num])
                    break
            if tail == eval(angle[2]):
                #print(head, tail, angle[2:], dihedral,3)
                t1 = eval(angle[3])
                t2 = eval(angle[4])
                if atoms[t1-1].symbol == dihedral[3] and atoms[t2-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum - num2, realxmaxnum - num,\
                                    t1 - num2, t2 - num2])
                    break
                elif atoms[t1-1].symbol == dihedral[2] and atoms[t2-1].symbol == dihedral[1]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum - num2, realxmaxnum - num,\
                                    t1 - num, t2 - num])
                    break
            if tail == eval(angle[4]):
                #print(head, tail, angle[2:], dihedral,4)
                t1 = eval(angle[2])
                t2 = eval(angle[3])
                if atoms[t2-1].symbol == dihedral[3] and atoms[t1-1].symbol == dihedral[4]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum - num2, realxmaxnum - num,\
                                    t2 - num, t1 - num])
                    break
                elif atoms[t2-1].symbol == dihedral[2] and atoms[t1-1].symbol == dihedral[1]:
                    #print(head, tail, angle[2:])
                    results.append([dihedral[0], realxminnum - num2, realxmaxnum - num,\
                                    t2 - num, t1 - num])
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
    if head - 1 in extra_angles[0]:
        head = head - 1
    elif head + 1 - num in extra_angles[0]:
        head = head + 1 - num
    if tail - 1 in extra_angles[0]:
        tail = tail - 1
    elif tail + 1 - num in extra_angles[0]:
        tail = tail + 1 - num
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
                    #print(angle[1:],angle2[1:],t, improper, [head,tail])
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
                    #print(angle[1:],angle2[1:],t, improper, [head,tail])
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

def bichain(file, degree, loc='./', show_bondtype=False, savename=''):
    """
    The method to create a single molecule for a finite chain based on a 2-monomer chain, specialised for polymers.
    bichain(file, degree, loc='./', savename='')
    file: Your molecule system name on your .data file.
    degree: Degree of polymerisation.
    loc: File Location. The default is your current location.
    savename: The name of the created file. The default is with '_XXx' (XX means the degree of polymerisation)
    Example:
        Input:
            from MCPoly.lmpset import chain
            bichain('BioPolymer2', 16)
        Output in Poly1_16x.data:
            1314 atoms
            1345 bonds
            2468 angles
            3382 dihedrals
            753 impropers
            
            9 atom types
            12 bond types
            21 angle types
            31 dihedral types
            3 improper types
            
            ...
            82 1 83 1
            
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
        a6 = re.search(r' xlo xhi', line)
        a7 = re.search(r' ylo yhi', line)
        a8 = re.search(r' zlo zhi', line)
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            tipnum = eval(b1.group(0))
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            tipnum2 = eval(b2.group(0))
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            tipnum3 = eval(b3.group(0))
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            tipnum4 = eval(b4.group(0))
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            tipnum5 = eval(b5.group(0))
        elif a6:
            b6 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            tipgapx = eval(b6[1]) - eval(b6[0])
            xlo = eval(b6[0])
            xhi = eval(b6[1])
        elif a7:
            b7 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapy = eval(b7[1]) - eval(b7[0])
            ylo = eval(b7[0])
        elif a8:
            b8 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapz = eval(b8[1]) - eval(b8[0])
            zlo = eval(b8[0])
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
                if eval(l1[0]) == tipnum:
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
                    Hloc[1] = [eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])]
                if eval(l[i4][-3]) < xmin:
                    xmin = eval(l[i4][-3])
                    xminnum = eval(l[i4][0])
                    Hloc[0] = [eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])]
            for i4 in range(len(l)):
                if eval(l[i4][2]) in H:
                    continue
                if eval(l[i4][-3]) > realxmax:
                    realxmax = eval(l[i4][-3])
                    realxmaxnum = eval(l[i4][0])
                if eval(l[i4][-3]) < realxmin:
                    realxmin = eval(l[i4][-3])
                    realxminnum = eval(l[i4][0])
            tipgapx = Hloc[1][0] - Hloc[0][0]
            w1 = 0
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
                    w2 = 2
                    j = 1
        elif w2 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost2_d = [eval(l[i4][-2]), eval(l[i4][-1])]
                    elif eval(l[i4][2]) == xmaxnum:
                        lost2_u = [eval(l[i4][-2]), eval(l[i4][-1])]
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost2_d = [eval(l[i4][-2]), eval(l[i4][-1])]
                    elif eval(l[i4][3]) == xmaxnum:
                        lost2_u = [eval(l[i4][-2]), eval(l[i4][-1])]
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
                if eval(l1[0]) == tipnum3:
                    w3 = 2
                    j = 1
        elif w3 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    if eval(l[i4][3]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost3_d.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
                    if eval(l[i4][4]) == xmaxnum:
                        lost3_u.append([eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1])])
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
                if eval(l1[0]) == tipnum4:
                    w4 = 2
                    j = 1
        elif w4 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][3]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][4]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][5]) in [xminnum, xmaxnum]:
                    if eval(l[i4][5]) == xminnum:
                        lost4_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][5]) == xmaxnum:
                        lost4_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
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
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            for i4 in range(len(l)):
                if eval(l[i4][2]) in [xminnum, xmaxnum]:
                    if eval(l[i4][2]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][2]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][3]) in [xminnum, xmaxnum]:
                    if eval(l[i4][3]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][3]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][4]) in [xminnum, xmaxnum]:
                    if eval(l[i4][4]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][4]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                elif eval(l[i4][5]) in [xminnum, xmaxnum]:
                    if eval(l[i4][5]) == xminnum:
                        lost5_d.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    elif eval(l[i4][5]) == xmaxnum:
                        lost5_u.append([eval(l[i4][-4]), eval(l[i4][-3]),\
                                        eval(l[i4][-2]), eval(l[i4][-1])])
                    continue
                j = j + 1
            k = k + 1
            l = []
            t5 = j
            j = 1
            k = 0
            w4 = 0
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    try:
        med = open('{0}_Chain.data'.format(file), 'r')
    except:
        infchain(file, stat[2], loc)
        med = open('{0}_Chain.data'.format(file), 'r')
    for line in med:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
        a6 = re.search(r' xlo xhi', line)
        
        m1 = re.search(r'Masses', line)
        if a1:
            b1 = re.search(r'[0-9]+', line)
            num = eval(b1.group(0))
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            num2 = eval(b2.group(0))
        elif a3:
            b3 = re.search(r'[0-9]+', line)
            num3 = eval(b3.group(0))
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            num4 = eval(b4.group(0))
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            num5 = eval(b5.group(0))
        elif a6:
            b6 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            gapx = eval(b6[1]) - eval(b6[0])
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
            if eval(l1[-1])<1.1:
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
                if eval(l1[0]) == tipnum:
                    w1 = 2
        elif w1 == 2:
            w1 = 0
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
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
                if eval(l1[0]) == tipnum3:
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
                if eval(l1[0]) == tipnum4:
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
                if eval(l1[0]) == tipnum5:
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
    try:
        if savename == '':
            f = open('{0}_{1}x.data'.format(file, degree*2), 'x')
        else:
            f = open('{0}.data'.format(savename), 'x')
    except:
        if savename == '':
            f = open('{0}_{1}x.data'.format(file, degree*2), 'w')
        else:
            f = open('{0}.data'.format(savename), 'w')
    inp.close()
    k = 0
    w1 = 0
    w2 = 0
    w3 = 0
    w4 = 0
    w5 = 0
    gap = 0
    l = []
    all_x = []
    inp = open('{0}.data'.format(file), 'r')
    print('X', xminnum, xmaxnum)
    print('REAL', realxminnum, realxmaxnum)
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
            f.write('{0} atoms\n'.format((tipnum-2)*degree + 2))
            #allnum = b1.group(0)
        elif a2:
            b2 = re.search(r'[0-9]+', line)
            f.write('{0} bonds\n'.format(tipnum2*degree - (degree-1)))
            #allnum2 = b2.group(0)
        elif a3:
            f.write('{0} angles\n'.format(tipnum3*degree))
            #allnum3 = b3.group(0)
        elif a4:
            b4 = re.search(r'[0-9]+', line)
            f.write('{0} dihedrals\n'.format(tipnum4*degree + 9*(degree-1)))
            #allnum4 = b4.group(0)
        elif a5:
            b5 = re.search(r'[0-9]+', line)
            f.write('{0} impropers\n'.format(tipnum5*degree))
            #allnum5 = b5.group(0)
        elif a6:
            x1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            xd = eval(x1[1]) - eval(x1[0])
            f.write('{0:>10.5f} {1:>10.5f} xlo xhi\n'.format(0.0,\
                                                       tipgapx*(degree-1) + gapx))
        elif a7:
            y1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            yd = eval(y1[1]) - eval(y1[0])
            f.write('{0:>10.5f} {1:>10.5f} ylo yhi\n'.format(0.0, gapy))
        elif a8:
            z1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            zd = eval(z1[1]) - eval(z1[0])
            f.write('{0:>10.5f} {1:>10.5f} zlo zhi\n'.format(0.0, gapz))
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
                if eval(l1[0]) == tipnum:
                    w1 = 2
                    j = 1
        elif w1 == 2:
            distance = 0
            allnum = 0
            for i3 in range(degree):
                for i4 in range(eval(b1.group(0))):
                    if i3 != 0 and i3 != degree-1:
                        if i4 + 1 in [xminnum, xmaxnum]:
                            gap = gap + 1
                            j = j + 1
                            continue
                    elif i3 == 0:
                        if i4 + 1 == xmaxnum:
                            gap = gap + 1
                            j = j + 1
                            all_x.append(eval(l[i4][-3]))
                            continue
                    elif i3 == degree - 1:
                        if i4 + 1 == xminnum:
                            gap = gap + 1
                            j = j + 1
                            continue
                    #print(i3, j-gap, eval(l[i4][-3]), eval(l[i4][-2]), eval(l[i4][-1]))
                    #if realxmin< = realxmax:
                    f.write('{0} {1} '.format(j-gap, eval(l[i4][1])))
                    f.write('{0} {1:>7.4f} '.format(eval(l[i4][2]), eval(l[i4][3])))
                    f.write('{0:>10.5f} '.format(eval(l[i4][-3]) + distance - xlo - (xmin-xlo)))
                    f.write('{0:>10.5f} '.format(eval(l[i4][-2]) - ylo))
                    f.write('{0:>10.5f}\n'.format(eval(l[i4][-1]) - zlo))
                    #else:
                    #    f.write('{0} {1} {2} {3:.4f} {4:.5f} {5:.5f} {6:.5f}\n'.format(j-gap,\
                    #eval(l[i4][1]), eval(l[i4][2]), eval(l[i4][3]),
                    #eval(l[i4][-3]) + distance-xhi, 
                    #eval(l[i4][-2])-ylo, eval(l[i4][-1])-zlo))
                    j = j + 1
                    if i3 == 0:
                        all_x.append(eval(l[i4][-3]))
                if i3 != 0 and i3 != degree - 1:
                    distance = distance + tipgapx
                    allnum = allnum + num - 2
                else:
                    distance = distance + gapx
                    allnum = allnum + num - 1
            l = []
            j = 1
            w1 = 0
            f.write('\n')
        elif w2 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 4:
                l.append(l1)
                if eval(l1[0]) == tipnum2:
                    lx = l
                    w2 = 2
                    j = 1
        elif w2 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(eval(b2.group(0))):
                    if i3 != 0 and i3 != degree - 1:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) - 2)
                        elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) - 2)
                        elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                    elif i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                    elif i3 == degree - 1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                    if i3 != 0 and i3 != degree - 1:
                        if eval(l[i4][2]) > eval(l[i4][3]):
                            f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                        else:
                            f.write('{0} {1} {3} {2}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                        j = j + 1
                    else:
                        if eval(l[i4][2]) <= tipnum-1 and eval(l[i4][3]) <= tipnum-1:
                            if eval(l[i4][2]) > eval(l[i4][3]):
                                f.write('{0} {1} {2} {3}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                            else:
                                f.write('{0} {1} {3} {2}\n'.format(j, eval(l[i4][1]),\
                                                                   eval(l[i4][2]) + every, eval(l[i4][3]) + every))
                            j = j + 1
                if i3 != 0:
                    if realxminnum <= realxmaxnum:
                        if realxminnum > realxmaxnum - num - 1:
                            f.write('{0} {1} {2} {3}\n'.format(j, bonds[0], realxminnum + every,\
                                                               realxmaxnum + every - num - 1))
                        else:
                            f.write('{0} {1} {3} {2}\n'.format(j, bonds[0], realxminnum + every,\
                                                               realxmaxnum + every - num - 1))
                    else:
                        if i3 == degree - 1:
                            if realxminnum > realxmaxnum - num:
                                f.write('{0} {1} {2} {3}\n'.format(j, bonds[0], realxminnum + every,\
                                                                   realxmaxnum + every - num))
                            else:
                                f.write('{0} {1} {3} {2}\n'.format(j, bonds[0], realxminnum + every,\
                                                                   realxmaxnum + every - num))
                        elif i3 == 1:
                            if realxminnum > realxmaxnum - num:
                                f.write('{0} {1} {2} {3}\n'.format(j, bonds[0], realxminnum + every - 1,\
                                                                   realxmaxnum + every - num - 1))
                            else:
                                f.write('{0} {1} {3} {2}\n'.format(j, bonds[0], realxminnum + every - 1,\
                                                                   realxmaxnum + every - num - 1))
                        else:
                            if realxminnum - 1 > realxmaxnum - num:
                                f.write('{0} {1} {2} {3}\n'.format(j, bonds[0], realxminnum + every - 1,\
                                                                   realxmaxnum + every - num))
                            else:
                                f.write('{0} {1} {3} {2}\n'.format(j, bonds[0], realxminnum + every - 1,\
                                                                   realxmaxnum + every - num))
                    j = j + 1
                #l = lx
                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum-1
                else:
                    every = every + tipnum-2
                for i4 in range(eval(b2.group(0))):
                    if i3 != 0 and i3 != degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        elif eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)
                        elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)
                        elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                    elif i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)  
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                    elif i3 == degree-1:
                        if (([eval(l[i4][2]), eval(l[i4][3])] == lost2_u
                                 or [eval(l[i4][2]), eval(l[i4][3])] == lost2_d) and
                            [eval(l[i4][2]), eval(l[i4][3])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                #if i3 == 1:
                    #print(l)
            k = k + 1
            all_bonds = l
            l = []
            j = 1
            w2 = 0
            k = 0
            f.write('\n')
        elif w3 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 5:
                l.append(l1)
                if eval(l1[0]) == tipnum3:
                    w3 = 2
                    j = 1
            #lx = l
        elif w3 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1) 
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                    elif i3 == degree-1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1) 
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3])-1)
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4])-1)           
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-2)  
                        elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-1)
                        if eval(l[i4][4]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-2)  
                        elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-1)
                    if eval(l[i4][2]) < eval(l[i4][4]):
                        f.write('{0} {1} {2} {3} {4}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every, eval(l[i4][4]) + every))
                    else:
                        f.write('{0} {1} {4} {3} {2}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][4]) + every, eval(l[i4][3]) + every, eval(l[i4][2]) + every))
                    j = j + 1
                if i3 == 1:
                    extra_angles = find_angle(file, head, tail, angles, all_bonds,\
                                              xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_angles:
                        rami_box = [*extra]
                        if extra[1] < 0 and extra[1] + tipnum > xmaxnum:
                            extra[1] = extra[1] - 1
                        if extra[2] < 0 and extra[2] + tipnum > xmaxnum:
                            extra[2] = extra[2] - 1
                        if extra[3] < 0 and extra[3] + tipnum > xmaxnum:
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
                            f.write('{0} {1} {2} {3} {4}\n'.format(j, extra[0], extra[1] + every,\
                                                                   extra[2] + every, extra[3] + every))
                        else:
                            f.write('{0} {1} {4} {3} {2}\n'.format(j, extra[0], extra[1] + every,\
                                                                   extra[2] + every, extra[3] + every))
                        extra[1] = rami_box[1]
                        extra[2] = rami_box[2]
                        extra[3] = rami_box[3]
                        j = j + 1
                elif i3 == degree - 1:
                    extra_angles = find_angle(file, head, tail, angles, all_bonds,\
                                              xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_angles:
                        rami_box = [*extra]
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] == min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] == min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] == min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        if extra[1] >= max([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] >= max([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] >= max([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[1] < extra[3]:
                            f.write('{0} {1} {2} {3} {4}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                   extra[2] + every + 1, extra[3] + every + 1))
                        else:
                            f.write('{0} {1} {4} {3} {2}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                   extra[2] + every + 1, extra[3] + every + 1))
                        extra[1] = rami_box[1]
                        extra[2] = rami_box[2]
                        extra[3] = rami_box[3]
                        j = j + 1
                elif i3 > 1:
                    extra_angles = find_angle(file, head, tail, angles, all_bonds,\
                                              xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_angles:
                        rami_box = [*extra]
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] >= min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] >= min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] >= min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        #if (extra[1] not in [min([realxminnum, realxmaxnum]), max([realxminnum, realxmaxnum])-1]
                        #        and extra[1] >= min([xminnum, xmaxnum])):
                        #    extra[1] = extra[1] - 1
                        #if (extra[2] not in [min([realxminnum, realxmaxnum]), max([realxminnum, realxmaxnum])-1]
                        #        and extra[2] >= min([xminnum, xmaxnum])):
                        #    extra[2] = extra[2] - 1
                        #if (extra[3] not in [min([realxminnum, realxmaxnum]), max([realxminnum, realxmaxnum])-1]
                        #        and extra[3] >= min([xminnum, xmaxnum])):
                        #    extra[3] = extra[3] - 1
                        if extra[1] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[3] = extra[3] - 1
                        #if (extra[1] not in [realxminnum, realxmaxnum]
                        #        and 0 < extra[1] < min([xminnum, xmaxnum])):
                        #    extra[1] = extra[1] + 1
                        #if (extra[2] not in [realxminnum, realxmaxnum]
                        #        and 0 < extra[2] < min([xminnum, xmaxnum])):
                        #    extra[2] = extra[2] + 1
                        #if (extra[3] not in [realxminnum, realxmaxnum]
                        #        and 0 < extra[3] < min([xminnum, xmaxnum])):
                        #    extra[3] = extra[3] + 1
                        if extra[1] < extra[3]:
                            f.write('{0} {1} {2} {3} {4}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                   extra[2] + every + 1, extra[3] + every + 1))
                        else:
                            f.write('{0} {1} {4} {3} {2}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                   extra[2] + every + 1, extra[3] + every + 1))
                        extra[1] = rami_box[1]
                        extra[2] = rami_box[2]
                        extra[3] = rami_box[3]
                        j = j + 1
                if i3 == 0 or i3 == degree-1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        elif eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1) 
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)  
                    elif i3 == degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        elif eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1) 
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)           
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_u
                                 or [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] in lost3_d) and
                            [eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]), eval(l[i4][4])])
                            continue
                        if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)
                        elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)
                        elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)
                        elif eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                #print(i3, l)
                #l = lx
            k = k + 1
            del all_bonds
            all_angles = l
            l = []
            j = 1
            k = 0
            w3 = 0
            f.write('\n')
        elif w4 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum4:
                    w4 = 2
                    j = 1
            #lx = l
        elif w4 == 2:
            every = 0
            for i3 in range(degree):
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)  
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        if eval(l[i4][5]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xmaxnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) - 1)
                    elif i3 == degree - 1:
                        if eval(l[i4][2]) == xminnum:
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)  
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)  
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        if eval(l[i4][5]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xminnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xminnum:
                            l[i4][5] = str(eval(l[i4][5]) - 1)
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            continue
                        elif eval(l[i4][2]) > max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-2)  
                        elif eval(l[i4][3]) > min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3])-1)
                        if eval(l[i4][4]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-2)  
                        elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4])-1)
                        if eval(l[i4][5]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-2)
                        elif eval(l[i4][5]) > min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5])-1)
                    f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every,\
                                                               eval(l[i4][4]) + every, eval(l[i4][5]) + every))
                    j = j + 1
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                    eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    elif i3 == degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xminnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] in lost4_u
                                or [eval(l[i4][2]), eval(l[i4][3]),\
                                    eval(l[i4][4]), eval(l[i4][5])] in lost4_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)
                        elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)
                        elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)
                        elif eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 2)
                        elif eval(l[i4][5]) >= min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                if i3 == 1:
                    extra_dihedrals = find_dihedral(file, head, tail, dihedrals, extra_angles,\
                                                    all_angles, xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_dihedrals:
                        if extra[1] < 0 and extra[1] + tipnum > xmaxnum:
                            extra[1] = extra[1] - 1
                        if extra[2] < 0 and extra[2] + tipnum > xmaxnum:
                            extra[2] = extra[2] - 1
                        if extra[3] < 0 and extra[3] + tipnum > xmaxnum:
                            extra[3] = extra[3] - 1
                        if extra[4] < 0 and extra[4] + tipnum > xmaxnum:
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
                        f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every,\
                                                                   extra[2] + every, extra[3] + every, extra[4] + every))
                        j = j + 1
                elif i3 == degree - 1:
                    extra_dihedrals = find_dihedral(file, head, tail, dihedrals, extra_angles,\
                                                    all_angles, xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_dihedrals:
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] <= 0 and extra[4] + tipnum > min([xminnum, xmaxnum]):
                            extra[4] = extra[4] - 1
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] == min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] == min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] == min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        if (extra[4] not in [realxminnum, realxmaxnum]
                                and extra[4] == min([xminnum, xmaxnum])):
                            extra[4] = extra[4] - 1
                        if extra[1] >= max([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] >= max([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] >= max([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] >= max([xminnum, xmaxnum]):
                            extra[4] = extra[4] - 1
                        f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                               extra[2] + every + 1, extra[3] + every + 1,\
                                                               extra[4] + every + 1))
                        j = j + 1
                elif i3 > 1:
                    extra_dihedrals = find_dihedral(file, head, tail, dihedrals, extra_angles,\
                                                    all_angles, xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for extra in extra_dihedrals:
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] <= 0 and extra[4] + tipnum > min([xminnum, xmaxnum]):
                            extra[4] = extra[4] - 1
                        #if (extra[1] not in [realxminnum, realxmaxnum]
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] >= min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] >= min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] >= min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        if (extra[4] not in [realxminnum, realxmaxnum]
                                and extra[4] >= min([xminnum, xmaxnum])):
                            extra[4] = extra[4] - 1
                        if extra[1] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[4] = extra[4] - 1
                        f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                   extra[2] + every + 1, extra[3] + every + 1,\
                                                                   extra[4] + every + 1))
                        j = j + 1
                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                #l = lx
            k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            f.write('\n')
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
        elif w5 == 2:
            lx = [*l]
            every = 0
            for i3 in range(degree):
                ignh = 0
                ignt = 0
                repeated = []
                for i4 in range(len(l)):
                    if i3 == 0:
                        #print(l)
                        if eval(l[i4][2]) == xmaxnum:
                            continue
                        elif eval(l[i4][2]) > xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            continue
                        elif eval(l[i4][3]) > xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        if eval(l[i4][4]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            continue
                        elif eval(l[i4][4]) > xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        if eval(l[i4][5]) == xmaxnum:
                            if eval(l[i4][2]) >= xmaxnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xmaxnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xmaxnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            continue
                        elif eval(l[i4][5]) > xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) - 1)
                        #print('1,end', l, [xminnum, xmaxnum])
                    elif i3 == degree - 1:
                        #print(l)
                        if str(xmaxnum) in l[i4][1:]:
                            ignt = ignt + 1
                        if eval(lx[i4][2]) == xminnum:
                            ignh = ignh + 1
                            continue
                        elif eval(l[i4][2]) > xminnum:
                            l[i4][2] = str(eval(l[i4][2]) - 1)
                        if eval(l[i4][3]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            ignh = ignh + 1
                            continue
                        elif eval(l[i4][3]) > xminnum:
                            l[i4][3] = str(eval(l[i4][3]) - 1)
                        if eval(l[i4][4]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            ignh = ignh + 1
                            continue
                        elif eval(l[i4][4]) > xminnum:
                            l[i4][4] = str(eval(l[i4][4])-1)
                        if eval(l[i4][5]) == xminnum:
                            if eval(l[i4][2]) >= xminnum:
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) >= xminnum:
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) >= xminnum:
                                l[i4][4] = str(eval(l[i4][4]) + 1)
                            ignh = ignh + 1
                            continue
                        elif eval(l[i4][5]) > xminnum:
                            l[i4][5] = str(eval(l[i4][5])-1)
                        #print('-1,end', l, [xminnum, xmaxnum])
                    else:
                        if eval(l[i4][2]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) == xminnum:
                                ignh = ignh + 1
                            elif eval(l[i4][2]) == xmaxnum:
                                ignt = ignt + 1
                            continue
                        if eval(l[i4][2]) > max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-2)  
                        elif eval(l[i4][2]) > min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2])-1)
                        if eval(l[i4][3]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
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
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) == xminnum:
                                ignh = ignh + 1
                            elif eval(l[i4][4]) == xmaxnum:
                                ignt = ignt + 1
                            continue
                        elif eval(l[i4][4]) > max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) - 2)  
                        elif eval(l[i4][4]) > min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) - 1)
                        if eval(l[i4][5]) in [xminnum, xmaxnum]:
                            if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 2)
                            elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                                l[i4][2] = str(eval(l[i4][2]) + 1)
                            if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 2)
                            elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                                l[i4][3] = str(eval(l[i4][3]) + 1)
                            if eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 2)
                            elif eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                                l[i4][4] = str(eval(l[i4][4]) + 1)
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
                                                               eval(l[i4][2]) + every, eval(l[i4][3]) + every,\
                                                               eval(l[i4][4]) + every, eval(l[i4][5]) + every))
                    j = j + 1
                for i4 in range(len(l)):
                    if i3 == 0:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                  or [eval(l[i4][2]), eval(l[i4][3]),\
                                      eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xmaxnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xmaxnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xmaxnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xmaxnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    elif i3 == degree - 1:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) >= xminnum:
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) >= xminnum:
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) >= xminnum:
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) >= xminnum:
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                    else:
                        if (([eval(l[i4][2]), eval(l[i4][3]),\
                              eval(l[i4][4]), eval(l[i4][5])] in lost5_u
                                 or [eval(l[i4][2]), eval(l[i4][3]),\
                                     eval(l[i4][4]), eval(l[i4][5])] in lost5_d) and
                            [eval(l[i4][2]), eval(l[i4][3]),\
                             eval(l[i4][4]), eval(l[i4][5])] not in repeated):
                            repeated.append([eval(l[i4][2]), eval(l[i4][3]),\
                                             eval(l[i4][4]), eval(l[i4][5])])
                            continue
                        if eval(l[i4][2]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 2)
                        elif eval(l[i4][2]) >= min([xminnum, xmaxnum]):
                            l[i4][2] = str(eval(l[i4][2]) + 1)
                        if eval(l[i4][3]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 2)
                        elif eval(l[i4][3]) >= min([xminnum, xmaxnum]):
                            l[i4][3] = str(eval(l[i4][3]) + 1)
                        if eval(l[i4][4]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 2)
                        elif eval(l[i4][4]) >= min([xminnum, xmaxnum]):
                            l[i4][4] = str(eval(l[i4][4]) + 1)
                        if eval(l[i4][5]) + 1 >= max([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 2)
                        elif eval(l[i4][5]) >= min([xminnum, xmaxnum]):
                            l[i4][5] = str(eval(l[i4][5]) + 1)
                #if i3 != 0:
                #    extra_impropers = find_improper(file, head, tail, impropers,\
                #                                    extra_angles, all_angles,\
                #                                    xminnum, xmaxnum, realxminnum, realxmaxnum)
                #    #print('ignh: ', ignh)
                #    #print('ignt: ', ignt, '\n')
                #    for z, extra in enumerate(extra_impropers):
                #        if z < len(extra_impropers) // 2:
                #            if z < len(extra_impropers) // 2 - ignh:
                #                continue
                #            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every,\
                #                                                       extra[2] + every, extra[3] + every, extra[4] + every))
                #        else:
                #            if z < len(extra_impropers) - ignt:
                #                continue
                #            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every,\
                #                                                       extra[2] + every, extra[3] + every, extra[4] + every))
                #        j = j + 1
                if i3 == 1:
                    print('ignh: ', ignh)
                    print('ignt: ', ignt)
                    extra_impropers = find_improper(file, head, tail, impropers,\
                                                     extra_angles, all_angles,\
                                                     xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for z,extra in enumerate(extra_impropers):
                        if extra[1] < 0 and extra[1] + tipnum > xmaxnum:
                            extra[1] = extra[1] - 1
                        if extra[2] < 0 and extra[2] + tipnum > xmaxnum:
                            extra[2] = extra[2] - 1
                        if extra[3] < 0 and extra[3] + tipnum > xmaxnum:
                            extra[3] = extra[3] - 1
                        if extra[4] < 0 and extra[4] + tipnum > xmaxnum:
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
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every,\
                                                                       extra[2] + every, extra[3] + every, extra[4] + every))
                        else:
                            if z < len(extra_impropers) - ignt:
                                continue
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every,\
                                                                       extra[2] + every, extra[3] + every, extra[4] + every))
                        j = j + 1
                elif i3 == degree - 1:
                    print('ignh: ', ignh)
                    print('ignt: ', ignt)
                    extra_impropers = find_improper(file, head, tail, impropers,\
                                                     extra_angles, all_angles,\
                                                     xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for z,extra in enumerate(extra_impropers):
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] <= 0 and extra[4] + tipnum > min([xminnum, xmaxnum]):
                            extra[4] = extra[4] - 1
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] == min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] == min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] == min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        if (extra[4] not in [realxminnum, realxmaxnum]
                                and extra[4] == min([xminnum, xmaxnum])):
                            extra[4] = extra[4] - 1
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
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                       extra[2] + every + 1, extra[3] + every + 1, extra[4] + every + 1))
                        else:
                            if z < len(extra_impropers) - ignt:
                                continue
                            print(z, len(extra_impropers), ignt)
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                       extra[2] + every + 1, extra[3] + every + 1, extra[4] + every + 1))
                        j = j + 1
                elif i3 > 1:
                    extra_impropers = find_improper(file, head, tail, impropers,\
                                                     extra_angles, all_angles,\
                                                     xminnum, xmaxnum, realxminnum, realxmaxnum)
                    for z,extra in enumerate(extra_impropers):
                        if extra[1] <= 0 and extra[1] + tipnum > min([xminnum, xmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] <= 0 and extra[2] + tipnum > min([xminnum, xmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] <= 0 and extra[3] + tipnum > min([xminnum, xmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] <= 0 and extra[4] + tipnum > min([xminnum, xmaxnum]):
                            extra[4] = extra[4] - 1
                        #if (extra[1] not in [realxminnum, realxmaxnum]
                        if (extra[1] not in [realxminnum, realxmaxnum]
                                and extra[1] >= min([xminnum, xmaxnum])):
                            extra[1] = extra[1] - 1
                        if (extra[2] not in [realxminnum, realxmaxnum]
                                and extra[2] >= min([xminnum, xmaxnum])):
                            extra[2] = extra[2] - 1
                        if (extra[3] not in [realxminnum, realxmaxnum]
                                and extra[3] >= min([xminnum, xmaxnum])):
                            extra[3] = extra[3] - 1
                        if (extra[4] not in [realxminnum, realxmaxnum]
                                and extra[4] >= min([xminnum, xmaxnum])):
                            extra[4] = extra[4] - 1
                        if extra[1] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[1] = extra[1] - 1
                        if extra[2] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[2] = extra[2] - 1
                        if extra[3] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[3] = extra[3] - 1
                        if extra[4] + 1 >= max([realxminnum, realxmaxnum]):
                            extra[4] = extra[4] - 1
                        if z < len(extra_impropers) // 2:
                            if z < len(extra_impropers) // 2 - ignh:
                                continue
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                       extra[2] + every + 1, extra[3] + every + 1, extra[4] + every + 1))
                        else:
                            if z < len(extra_impropers) - ignt:
                                continue
                            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, extra[0], extra[1] + every + 1,\
                                                                       extra[2] + every + 1, extra[3] + every + 1, extra[4] + every + 1))
                        j = j + 1
                if i3 == 0 or i3 == degree - 1:
                    every = every + tipnum - 1
                else:
                    every = every + tipnum - 2
                #lx = l
            k = k + 1
            l = []
            j = 1
            k = 0
            w4 = 0
            l = lx
        else:
            f.write(line)
    f.close()
    inp.close()
    med.close()
    os.chdir(opath)
    if show_bondtype == True:
        return eval(bonds[0])