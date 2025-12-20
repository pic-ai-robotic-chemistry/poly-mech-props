import os
import warnings
from ase.io import read
import random as r
import sys

mydir = os.path.dirname( __file__ )
statusdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(statusdir)
from coordsys import coordsys

def distance(dot1, dot2):
    delta_x = (dot1[0] - dot2[0]) ** 2
    delta_y = (dot1[1] - dot2[1]) ** 2
    delta_z = (dot1[2] - dot2[2]) ** 2
    return (delta_x + delta_y + delta_z) ** 0.5

def randombuild(file, num, gap=2.0, size=[100, 100, 100], bg='empty', bgloc='./', loc='./', saveloc='./', savename='', errornum=100):
    if savename == '':
        savename = file
    savename = savename + '_{0}.xyz'.format(num)
    opath = os.getcwd()
    w = open(saveloc+savename, 'w')
    try:
        atoms_zmx = coordsys(file, filefmt='xyz', coordtype='zmx', loc=loc, save=False)
    except:
        atoms_zmx = coordsys(file, filefmt='zmx', coordtype='zmx', loc=loc, save=False, warn=False)
    i = 0
    step = 0
    base = []
    if bg != 'empty':
        atoms_xyz0 = coordsys(bg, filefmt='xyz', coordtype='xyz', loc=bgloc, save=False, warn=False)
        base = atoms_xyz0
        w.write('{0}\n\n'.format(len(atoms_zmx)*num + len(base)))
        for coord in base:
            #print(coord)
            coord[1] = eval(coord[1])
            coord[2] = eval(coord[2])
            coord[3] = eval(coord[3])
    else:
        w.write('{0}\n\n'.format(len(atoms_zmx)*num))
    base_pre = []
    judge = True
    while 1:
        start_point = [r.random()*size[0], r.random()*size[1], r.random()*size[2]]
        roto = [r.random()*360-180, r.random()*360-180, r.random()*360-180]
        atoms_xyz = coordsys(file, filefmt='zmx', coordtype='xyz', loc=loc, rotate=roto, savename='mid')
        #print(atoms_xyz)
        for atom in atoms_xyz:
            a = atom[1] + start_point[0]
            b = atom[2] + start_point[1]
            c = atom[3] + start_point[2]
            dot1 = [a, b, c]
            #print(base)
            for coord in base:
                d = distance(dot1, [coord[1], coord[2], coord[3]])
                #print(dot1,[coord[1], coord[2], coord[3]])
                if d >= gap:
                    judge = True
                else:
                    judge = False
                    break
            if dot1[0] > size[0] or dot1[0] < 0:
                judge = False
            if dot1[1] > size[1] or dot1[1] < 0:
                judge = False
            if dot1[2] > size[2] or dot1[2] < 0:
                judge = False
            if judge == False:
                base_pre = []
                break
            else:
                base_pre.append([atom[0], *dot1])
        if judge == True:
            i = i + 1
            base = base + base_pre
            base_pre = []
        else:
            base_pre = []
            judge = True
        step = step + 1
        if i == num:
            break
        if step > num + errornum:
            os.system('rm mid.xyz')
            raise RuntimeError('Time Runs out.')
    for atom in base:
        w.write('{0:<2} {1:>9.5f} {2:>9.5f} {3:>9.5f}\n'.format(atom[0], atom[1], atom[2], atom[3]))
    w.close()
    os.system('rm mid.xyz')
    os.chdir(opath)