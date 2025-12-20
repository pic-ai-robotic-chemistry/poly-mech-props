import re
import os
import sys

mydir = os.path.dirname( __file__ )
lmpdir = os.path.join(mydir, '..', 'lmpset')
sys.path.append(lmpdir)
from infchain import infchain

def order(file, orderlist, loc='./', savename=''):
    """
    The method to recreate a SINGLE FINITE molecule in a different atom order.
    order(file, bondtype, degree, loc='./', savename='')
    file: Your molecule system name on your .data file.
    orderlist: The new order of all atoms.
    loc: File Location. The default is your current location.
    savename: The name of the created file. The default is with '_od'
    Example:
        Input:
            from MCPoly.lmpset import chain
            chain('BioPolymer3', 11, 16)
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
    tipnum = len(orderlist)
    orderdict = {}
    for i,unit in enumerate(orderlist):
        orderdict[unit+1] = i + 1
    print(orderdict)
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
    real = [0, 0, 0]
    inp = open('{0}.data'.format(file), 'r')
    for line in inp:
        a1 = re.search(r' atoms', line)
        a2 = re.search(r' bonds', line)
        a3 = re.search(r' angles', line)
        a4 = re.search(r' dihedrals', line)
        a5 = re.search(r' impropers', line)
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
    try:
        if savename == '':
            f = open('{0}_od.data'.format(file), 'x')
        else:
            f = open('{0}.data'.format(savename), 'x')
    except:
        if savename == '':
            f = open('{0}_od.data'.format(file), 'w')
        else:
            f = open('{0}.data'.format(savename), 'w')
    inp.close()
    inp = open('{0}.data'.format(file), 'r')
    l = []
    all_x = []
    for line in inp:
        c1 = re.search('Atoms', line)
        c2 = re.search('Bonds', line)
        c3 = re.search('Angles', line)
        c4 = re.search('Dihedrals', line)
        c5 = re.search('Impropers', line)
        if c1:
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
            for i3, elenum in enumerate(orderlist):
                f.write('{0} {1} '.format(j, eval(l[elenum][1])))
                f.write('{0} {1:>7.4f} '.format(eval(l[elenum][2]), eval(l[elenum][3])))
                f.write('{0:>10.5f} '.format(eval(l[elenum][-3])))
                f.write('{0:>10.5f} '.format(eval(l[elenum][-2])))
                f.write('{0:>10.5f}\n'.format(eval(l[elenum][-1])))
                j = j + 1
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
            for tie in l:
                tie[0] = eval(tie[0])
                tie[1] = eval(tie[1])
                tie[-2] = orderdict[eval(tie[-2])]
                tie[-1] = orderdict[eval(tie[-1])]
            re_l = []
            j2 = 1
            for tie in l:
                if tie[2] < tie[3]:
                    x = tie[2]
                    tie[2] = tie[3]
                    tie[3] = x
            for i in range(tipnum):
                inter_l = []
                for tie in l:
                    if tie[-1] == i + 1:
                        inter_l.append(tie)
                if len(inter_l) > 1:
                    for j0 in range(tipnum):
                        for tie2 in inter_l:
                            if tie2[-2] == j0 + 1:
                                re_l.append([j2, *tie2[1:]])
                                j2 = j2 + 1
                elif len(inter_l) == 1:
                    re_l.append([j2, *inter_l[0][1:]])
                    j2 = j2 + 1
            for i4 in range(eval(b2.group(0))):
                f.write('{0} {1} {2} {3}\n'.format(j, re_l[i4][1],\
                                                   re_l[i4][2], re_l[i4][3]))
                j = j + 1
            l = []
            j = 1
            w2 = 0
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
            for tie in l:
                tie[-4] = eval(tie[-4])
                tie[-3] = orderdict[eval(tie[-3])]
                tie[-2] = orderdict[eval(tie[-2])]
                tie[-1] = orderdict[eval(tie[-1])]
                tie[0] = eval(tie[0])
            re_l = []
            j2 = 1
            for tie in l:
                if tie[2] > tie[4]:
                    x = tie[2]
                    tie[2] = tie[4]
                    tie[4] = x
            for i in range(tipnum):
                inter_l = []
                for tie in l:
                    if tie[-3] == i + 1:
                        inter_l.append(tie)
                if len(inter_l) > 1:
                    for j0 in range(tipnum):
                        inter_l2 = []
                        for tie2 in inter_l:
                            if tie2[-2] == j0 + 1:
                                inter_l2.append([j2, *tie2[1:]])
                        if len(inter_l2) > 1:
                            for k in range(tipnum):
                                for tie3 in inter_l2:
                                    if tie3[-2] == k + 1:
                                        re_l.append([j2, *tie3[1:]])
                                        j2 = j2 + 1
                        elif len(inter_l2) == 1:
                            re_l.append([j2, *inter_l2[0][1:]])
                            j2 = j2 + 1
                elif len(inter_l) == 1:
                    re_l.append([j2, *inter_l[0][1:]])
                    j2 = j2 + 1
            for i4 in range(eval(b3.group(0))):
                f.write('{0} {1} {2} {3} {4}\n'.format(j, re_l[i4][1], re_l[i4][2],\
                                                       re_l[i4][3], re_l[i4][4]))
                j = j + 1
            l = []
            j = 1
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
            for tie in l:
                tie[-4] = orderdict[eval(tie[-4])]
                tie[-3] = orderdict[eval(tie[-3])]
                tie[-2] = orderdict[eval(tie[-2])]
                tie[-1] = orderdict[eval(tie[-1])]
                tie[0] = eval(tie[0])
                tie[1] = eval(tie[1])
            re_l = []
            for i in range(tipnum):
                inter_l = []
                for tie in l:
                    if tie[-1] == i + 1:
                        re_l.append(tie)
            for tie in re_l:
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, *tie[1:]))
                j = j + 1
        #    every = 0
        #    for i3 in range(degree):
        #        repeated = []
        #        for i4 in range(len(l)):
        #            f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, eval(l[i4][1]),\
        #                                                       eval(l[i4][2]) + every, eval(l[i4][3]) + every,\
        #                                                       eval(l[i4][4]) + every, eval(l[i4][5]) + every))
        #            j = j + 1
        #    k = k + 1
            l = []
        #    j = 1
        #    k = 0
            w4 = 0
            f.write('\n')
        elif w5 == 1:
            l1 = re.findall(r'[0-9]+', line)
            if len(l1) == 6:
                l.append(l1)
                if eval(l1[0]) == tipnum5:
                    w5 = 2
                    j = 1
            #lx = l
        elif w5 == 2:
            for tie in l:
                tie[-4] = orderdict[eval(tie[-4])]
                tie[-3] = orderdict[eval(tie[-3])]
                tie[-2] = orderdict[eval(tie[-2])]
                tie[-1] = orderdict[eval(tie[-1])]
                tie[0] = eval(tie[0])
                tie[1] = eval(tie[1])
            re_l = []
            for i in range(tipnum):
                inter_l = []
                for tie in l:
                    if tie[-1] == i + 1:
                        re_l.append(tie)
            for tie in re_l:
                f.write('{0} {1} {2} {3} {4} {5}\n'.format(j, *tie[1:]))
                j = j + 1
        #    every = 0
        #    k = k + 1
            l = []
        #    j = 1
        #    k = 0
            w5 = 0
        #    l = lx
        #    f.write('\n')
        else:
            f.write(line)
    f.close()
    inp.close()
    os.chdir(opath)