import re
import os

def DumpToXyz(dumpfile, datafile='', loc='./', savefile=''):
    """
    The method to create a XYZ trajectory file based on LAMMPS DUMP file.
    DumpToXyz(dumpfile, datafile='', loc='./', savefile='')
    dumpfile: Your LAMMPS DUMP name on your .dump file.
    datafile: Your LAMMPS Data name on your .data file.
    loc: File Location. The default is your current location.
    datafile: Your name of your saved .XYZ trajectory file.
    """
    path = os.getcwd()
    os.chdir(loc)
    if datafile == '':
        datafile = dumpfile
    data = open('{0}.data'.format(datafile), 'r')
    if savefile == '':
        savefile = datafile
    elementsort = []
    multiplier = []
    l = []
    n1 = 0
    n2 = 0
    for line in data:
        m1 = re.search('Masses', line)
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
            mass = eval(l1[-1])
            if abs(mass - 1) < 0.5:
                elementsort.append('H')
            elif abs(mass - 2) < 0.5:
                elementsort.append('D')
            elif abs(mass - 3) < 0.5:
                elementsort.append('T')
            elif abs(mass - 4) < 0.5:
                elementsort.append('He')
            elif abs(mass - 7) < 0.5:
                elementsort.append('Li')
            elif abs(mass - 9) < 0.5:
                elementsort.append('Be')
            elif abs(mass - 10.8) < 0.5:
                elementsort.append('B')
            elif abs(mass - 12) < 0.5:
                elementsort.append('C')
            elif abs(mass - 14) < 0.5:
                elementsort.append('N')
            elif abs(mass - 16) < 0.5:
                elementsort.append('O')
            elif abs(mass - 19) < 0.5:
                elementsort.append('F')
            elif abs(mass - 22.9) < 0.5:
                elementsort.append('Na')
            elif abs(mass - 24.3) < 0.5:
                elementsort.append('Mg')
            elif abs(mass - 27) < 0.5:
                elementsort.append('Al')
            elif abs(mass - 28.1) < 0.5:
                elementsort.append('Si')
            elif abs(mass - 31) < 0.5:
                elementsort.append('P')
            elif abs(mass - 32) < 0.5:
                elementsort.append('S')
            elif abs(mass - 35.45) < 0.5:
                elementsort.append('Cl')
            elif abs(mass - 80) < 0.5:
                elementsort.append('Br')
            elif abs(mass - 127) < 0.5:
                elementsort.append('I')
            elif abs(mass - 39.1) < 0.5:
                elementsort.append('K')
            elif abs(mass - 40.1) < 0.5:
                elementsort.append('Ca')
            elif abs(mass - 47.86) < 0.2:
                elementsort.append('Ti')
            elif abs(mass - 54.93) < 0.2:
                elementsort.append('Mn')
            elif abs(mass - 55.84) < 0.2:
                elementsort.append('Fe')
            elif abs(mass - 58.93) < 0.12:
                elementsort.append('Co')
            elif abs(mass - 58.69) < 0.12:
                elementsort.append('Ni')
            elif abs(mass - 63.55) < 0.2:
                elementsort.append('Cu')
            elif abs(mass - 65.4) < 0.2:
                elementsort.append('Zn')
            elif abs(mass - 101.07) < 0.2:
                elementsort.append('Ru')
            elif abs(mass - 102.91) < 0.2:
                elementsort.append('Rh')
            elif abs(mass - 106.42) < 0.2:
                elementsort.append('Pd')
            elif abs(mass - 107.87) < 0.2:
                elementsort.append('Ag')
            elif abs(mass - 112.41) < 0.2:
                elementsort.append('Cd')
            elif abs(mass - 192.22) < 0.2:
                elementsort.append('Ir')
            elif abs(mass - 195.08) < 0.2:
                elementsort.append('Pt')
            elif abs(mass - 196.97) < 0.2:
                elementsort.append('Au')
            elif abs(mass - 200.59) < 0.2:
                elementsort.append('Hg')
        if m1:
            n1 = 1
    data.close()
    print(elementsort)
    f = open('{0}_dump.xyz'.format(savefile), 'w')
    dump = open('{0}.dump'.format(dumpfile), 'r')
    for line in dump:
        m1 = re.search('ITEM: ATOMS id type xs ys zs', line)
        m2 = re.search('ITEM: BOX BOUNDS pp pp pp', line)
        #m3 = re.search('ITEM: NUMBER OF ATOMS', line)
        if n2 == 1:
            x2 = re.search(r'.+[0-9] ', line)
            x3 = re.search(r' .+[0-9]', line)
            if x2:
                multiplier.append(eval(x3.group(0)[1:]) - eval(x2.group(0)[:-1]))
            else:
                n2 = 0
        if n1 == 1:
            l1 = re.findall(r'\-?[0-9]+\.?[0-9]*', line)
            if len(l1) != 5:
                length = len(l)
                n1 = 0
                f.write('{0}\n'.format(length))
                f.write('Dump Result\n')
                for figures in l:
                    f.write('  {0}       {1:>10.6f}       '.format(*figures[0:2]))
                    f.write('{0:>10.6f}       {1:>10.6f}\n'.format(*figures[2:]))
                l = []
            else:
                num1 = eval(l1[1]) - 1
                num2 = eval(l1[2])
                num3 = eval(l1[3])
                num4 = eval(l1[4])
                l.append([elementsort[num1], num2 * multiplier[0],\
                             num3 * multiplier[1], num4 * multiplier[2]])
        if m1:
            n1 = 1
        if m2:
            n2 = 1
    os.chdir(path)