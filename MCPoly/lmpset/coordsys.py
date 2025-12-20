import os
import warnings
import re
import math as m

def vector(atoms, atomnum1, atomnum2):
    x = eval(atoms[atomnum1][0]) - eval(atoms[atomnum2][0])
    y = eval(atoms[atomnum1][1]) - eval(atoms[atomnum2][1])
    z = eval(atoms[atomnum1][2]) - eval(atoms[atomnum2][2])
    return (x, y, z)

def dot_product(vector1, vector2):
    a = vector1[0] * vector2[0]
    b = vector1[1] * vector2[1]
    c = vector1[2] * vector2[2]
    return a + b + c

def cross_product(vector1, vector2, value=False):
    a = vector1[1] * vector2[2] - vector1[2] * vector2[1]
    b = vector1[2] * vector2[0] - vector1[0] * vector2[2]
    c = vector1[0] * vector2[1] - vector1[1] * vector2[0]
    if value == False:
        return (a, b, c)

def judge(vector1, vector2, vector0):
    test = cross_product(vector1, vector2)
    if test[0] * vector0[0] > 0.01:
        return 1
    elif test[0] * vector0[0] < -0.01:
        return -1
    else:
        if test[1] * vector0[1] > 0.01:
            return 1
        elif test[1] * vector0[1] < -0.01:
            return -1
        else:
            if test[2] * vector0[2] > 0.01:
                return 1
            elif test[2] * vector0[2] < -0.01:
                return -1

def rotate_x(x0, y0, z0, theta):
    x = x0
    y = y0 * m.cos(theta * m.pi / 180) + z0 * m.sin(theta * m.pi / 180)
    z = -y0 * m.sin(theta * m.pi / 180) + z0 * m.cos(theta * m.pi / 180)
    return [x, y, z]

def rotate_y(x0, y0, z0, theta):
    x = x0 * m.cos(theta * m.pi / 180) + z0 * m.sin(theta * m.pi / 180)
    y = y0
    z = -x0 * m.sin(theta * m.pi / 180) + z0 * m.cos(theta * m.pi / 180)
    return [x, y, z]

def rotate_z(x0, y0, z0, theta):
    x = x0 * m.cos(theta * m.pi / 180) + y0 * m.sin(theta * m.pi / 180)
    y = -x0 * m.sin(theta * m.pi / 180) + y0 * m.cos(theta * m.pi / 180)
    z = z0
    return [x, y, z]

def cartestance(atoms, atomnum1, atomnum2):
    delta_x = (eval(atoms[atomnum1][0]) - eval(atoms[atomnum2][0])) ** 2
    delta_y = (eval(atoms[atomnum1][1]) - eval(atoms[atomnum2][1])) ** 2
    delta_z = (eval(atoms[atomnum1][2]) - eval(atoms[atomnum2][2])) ** 2
    d = (delta_x + delta_y + delta_z) ** 0.5
    return d

def angle(atoms, atomnum1, atomnum2, atomnum3):
    d1 = cartestance(atoms, atomnum1, atomnum2)
    d2 = cartestance(atoms, atomnum2, atomnum3)
    vector1 = vector(atoms, atomnum2, atomnum1)
    vector2 = vector(atoms, atomnum2, atomnum3)
    prod = dot_product(vector1, vector2)
    cos_theta = prod / (d1 * d2)
    theta = m.acos(cos_theta) / m.pi * 180
    return theta

def dihedral(atoms, atomnum1, atomnum2, atomnum3, atomnum4):
    vector1 = vector(atoms, atomnum2, atomnum1)
    vector2 = vector(atoms, atomnum2, atomnum3)
    vector3 = vector(atoms, atomnum4, atomnum3)
    vertivector1 = cross_product(vector2, vector1)
    vertivector2 = cross_product(vector2, vector3)
    state = judge(vertivector1, vertivector2, vector2)
    d1 = vertivector1[0] ** 2 + vertivector1[1] ** 2 + vertivector1[2] ** 2
    d1 = d1 ** 0.5
    d2 = vertivector2[0] ** 2 + vertivector2[1] ** 2 + vertivector2[2] ** 2
    d2 = d2 ** 0.5
    prod = dot_product(vertivector1, vertivector2)
    cos_theta = prod / (d1 * d2)
    theta = m.acos(cos_theta) / m.pi * 180
    return state * (180 - theta)

def coordsys(file, filefmt='xyz', coordtype='zmx', loc='./', save=True, savename='', saveloc='./', rotate=[0, 0, 0], warn=True):
    theta_x = rotate[0]
    theta_y = rotate[1]
    theta_z = rotate[2]
    if savename == '':
        savename = file
    opath = os.getcwd()
    if filefmt == 'xyz':
        file = file + '.xyz'
    elif filefmt == 'zmx':
        file = file + '.zmx'
    f = open(loc+file, 'r')
    if save == True:
        if coordtype == 'xyz':
            try:
                w = open(saveloc+savename+'.xyz', 'x')
            except:
                w = open(saveloc+savename+'.xyz', 'w')
        elif coordtype == 'zmx':
            try:
                w = open(saveloc+savename+'.zmx', 'x')
            except:
                w = open(saveloc+savename+'.zmx', 'w')
    else:
        pass
    atoms = []
    alls = []
    eles = []
    newatoms = []
    newalls = []
    for i,line in enumerate(f):
        if i < 2:
            if save == True:
                w.write(line)
        elif i >= 2:
            ele = re.search(r'[A-Z][a-z]?', line)
            count = re.findall(r'\-?[0-9]+\.[0-9]+', line)
            for x in count:
                x = eval(x)
            atoms.append(count)
            eles.append(ele.group(0))
            alls.append([ele.group(0), *count])
            if filefmt == coordtype and save == True:
                w.write(line)
    if filefmt == coordtype and save == True:
        w.close()
    f.close()
    file = file + '.' + filefmt
    if filefmt == coordtype:
        if warn == True:
            print("You've just copied the file in the same format!")
        os.chdir(opath)
        return alls
    if filefmt == 'xyz' and coordtype == 'zmx':
        for i,atom in enumerate(atoms):
            if i == 0:
                newatoms.append([])
                newalls.append([eles[i]])
                if save == True:
                    w.write('{0}\n'.format(eles[i]))
            elif i == 1 :
                d1 = cartestance(atoms, i, 0)
                newatoms.append([d1])
                newalls.append([eles[i], d1])
                if save == True:
                    w.write('{0:<2} {1:>9.5f}\n'.format(eles[i], d1))
            elif i == 2 :
                d1 = cartestance(atoms, i, 0)
                angle1 = angle(atoms, i, 0, 1)
                newatoms.append([d1, angle1])
                newalls.append([eles[i], d1, angle1])
                if save == True:
                    w.write('{0:<2} {1:>9.5f} {2:>9.3f}\n'.format(eles[i], d1, angle1))
            else:
                d1 = cartestance(atoms, i, 0)
                angle1 = angle(atoms, i, 0, 1)
                dihedral1 = dihedral(atoms, i, 0, 1, 2)
                newatoms.append([d1, angle1, dihedral1])
                newalls.append([eles[i], d1, angle1, dihedral1])
                if save == True:
                    w.write('{0:<2} {1:>9.5f} {2:>9.3f} {3:>9.3f}\n'.format(eles[i], d1, angle1, dihedral1))
        if save == True:
            w.close()
        os.chdir(opath)
        return newalls
    elif filefmt == 'zmx' and coordtype == 'xyz':
        for i,atom in enumerate(atoms):
            if i == 0:
                newatoms.append([0.00000, 0.00000, 0.00000])
                newalls.append([eles[i], 0.00000, 0.00000, 0.00000])
                if save == True:
                    w.write('{0:<2} {1:>9.5f} {2:>9.5f} {3:>9.5f}\n'.format(eles[i], 0.00000, 0.00000, 0.00000))
            elif i == 1:
                a = eval(atoms[i][0])
                b = 0
                c = 0
                coord = rotate_x(a, b, c, theta_x)
                coord = rotate_y(*coord, theta_y)
                coord = rotate_z(*coord, theta_z)
                newatoms.append([coord])
                newalls.append([eles[i], *coord])
                if save == True:
                     w.write('{0:<2} {1:>9.5f} {2:>9.5f} {3:>9.5f}\n'.format(eles[i], *coord))
            elif i == 2:
                a = eval(atoms[i][0]) * m.cos(eval(atoms[i][1]) * m.pi/180)
                b = eval(atoms[i][0]) * m.sin(eval(atoms[i][1]) * m.pi/180)
                c = 0
                coord = rotate_x(a, b, c, theta_x)
                coord = rotate_y(*coord, theta_y)
                coord = rotate_z(*coord, theta_z)
                newatoms.append(coord)
                newalls.append([eles[i], *coord])
                if save == True:
                    w.write('{0:<2} {1:>9.5f} {2:>9.5f} {3:>9.5f}\n'.format(eles[i], *coord))
            else:
                a = eval(atoms[i][0]) * m.cos(eval(atoms[i][1]) * m.pi/180)
                b = eval(atoms[i][0]) * m.sin(eval(atoms[i][1]) * m.pi/180) * m.cos(eval(atoms[i][2]) * m.pi/180)
                c = -eval(atoms[i][0]) * m.sin(eval(atoms[i][1]) * m.pi/180) * m.sin(eval(atoms[i][2]) * m.pi/180)
                coord = rotate_x(a, b, c, theta_x)
                coord = rotate_y(*coord, theta_y)
                coord = rotate_z(*coord, theta_z)
                newatoms.append(coord)
                newalls.append([eles[i], *coord])
                if save == True:
                    w.write('{0:<2} {1:>9.5f} {2:>9.5f} {3:>9.5f}\n'.format(eles[i], *coord))
        if save == True:
            w.close()
        os.chdir(opath)
        return newalls