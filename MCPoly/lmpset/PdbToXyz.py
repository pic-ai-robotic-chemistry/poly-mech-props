import re

def PdbToXyz(file, loc='./', saveloc='./', savename=''):
    """
    The method to convert .pdb file into .xyz file.
    PdbToXyz(file, loc='./', saveloc='./', savename='')
    file: Your .pdb file name.
    loc: File Location. The default is your current location.
    saveloc: Input File Save Location. The default is your current location.
    savename: The name of your saved .xyz file. Default is the current name of .pdb file.
    Example 1:
        Input:
            from MCPoly.lmpset import PdbToXyz
            PdbToXyz('Mole1')
        Output in Mole1.xyz:
            336
            #Powered by MCPoly
            O         45.43700        9.79500       13.77600
            C         45.99400        9.03600       12.98700
            C         45.33800        8.04800       12.00700
            O         47.33200        9.06400       12.77100
    """
    f = open(loc+file+'.pdb', 'r')
    if savename == '':
        savename = file
    try:
        w = open(saveloc+savename+'.xyz', 'x')
    except:
        w = open(saveloc+savename+'.xyz', 'w')
    for line in f:
        a = re.search('ATOM', line)
        if a:
            i = re.search(r' [0-9]+ ', line)
        z = re.search('END', line)
        if z:
            w.write(i.group(0)[1:-1])
            w.write('\n#Powered by MCPoly\n')
            break
    f.close()
    f = open(loc+file+'.pdb', 'r')
    for line in f:
        a = re.search('ATOM', line)
        if a:
            ele = re.search(' [A-Z][a-z]* ', line)
            b = re.findall(r'\-?[0-9]+\.[0-9]+', line)
            w.write('  {0}       {1:>10.6f}       '.format(ele.group(0)[1:-1], eval(b[0])))
            w.write('{0:>10.6f}       {1:>10.6f}\n'.format(eval(b[1]), eval(b[2])))
        z = re.search('END', line)
        if z:
            w.write('\n')
            break
    w.close()