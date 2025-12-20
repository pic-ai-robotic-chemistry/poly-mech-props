import re

def PdbToXyz_Multi(file, steps=[], loc='./', saveloc='./', savename=''):
    """
    The method to convert .pdb trajectory file into .xyz single file and trajectory file.
    PdbToXyz_Multi(file, steps=[], loc='./', saveloc='./', savename='')
    file: Your .pdb file name.
    steps: The step number
    loc: File Location. The default is your current location.
    saveloc: Input File Save Location. The default is your current location.
    savename: The name of your saved .xyz file. Default is the current name of .pdb file.
    Example 1:
        Input:
            from MCPoly.lmpset import PdbToXyz_Multi
            PdbToXyz_Multi('MoleTraj1')
        Output in MoleTraj1.xyz:
            336
            #Powered by MCPoly
            O         45.43700        9.79500       13.77600
            C         45.99400        9.03600       12.98700
            C         45.33800        8.04800       12.00700
            O         47.33200        9.06400       12.77100
            ......
            With all structure in PDB file.
            
    Example 2:
        Input:
            from MCPoly.lmpset import PdbToXyz_Multi
            PdbToXyz_Multi('MoleTraj1',savename='MoleTraj2', steps=[530])
        Output in MoleTraj2.xyz:
            336         The Structure of Step 530
            #Powered by MCPoly
            O        45.353000         9.899000        13.941000
            C        45.918000         8.891000        13.494000
            C        45.215000         7.741000        12.699000
            O        47.237000         8.549000        13.612000

    Example 3:
        Input:
            from MCPoly.lmpset import PdbToXyz_Multi
            PdbToXyz_Multi('MoleTraj1',savename='MoleTraj2', steps=[-5,-4,-3,-2,-1])
        Output in MoleTraj2.xyz:
            336         The Structure of Step 530
            #Powered by MCPoly
            O        45.444000         9.785000        13.789000
            C        45.995000         9.028000        13.006000
            C        45.365000         8.054000        12.037000
            O        47.326000         9.041000        12.760000
            ...
            336
            #Powered by MCPoly
            O        45.442000         9.790000        13.784000
            C        45.995000         9.029000        12.998000
            C        45.355000         8.053000        12.026000
            O        47.328000         9.051000        12.764000
            ......
            Trajectory of the last five steps.
    """
    step = 1
    f = open(loc+file+'.pdb', 'r')
    if savename == '':
        savename = file
    try:
        w = open(saveloc+savename+'.xyz', 'x')
    except:
        w = open(saveloc+savename+'.xyz', 'w')
    stepnum = 0
    for line in f:
        a = re.search('ATOM', line)
        if a:
            i = re.search(r' [0-9]+ ', line)
        z = re.search('END', line)
        if z and stepnum == 0:
            if steps == [] or 1 in steps:
                w.write(i.group(0)[1:-1])
                w.write('\n#Powered by MCPoly\n')
            stepnum = 1
        elif z and stepnum != 0:
            stepnum = stepnum + 1
    f.close()
    f = open(loc+file+'.pdb', 'r')
    for line in f:
        a = re.search('ATOM', line)
        if a:
            ele = re.search(r' [A-Z][a-z]* ', line)
            b = re.findall(r'\-?[0-9]+\.[0-9]+', line)
            if (steps == []) or (step in steps) or (step-stepnum-1 in steps):
                w.write('  {0}       {1:>10.6f}       '.format(ele.group(0)[1:-1], eval(b[0])))
                w.write('{0:>10.6f}       {1:>10.6f}\n'.format(eval(b[1]), eval(b[2])))
        z = re.search('END', line)
        if z and step+1 <= stepnum:
            if (steps == []) or (step+1 in steps) or (step-stepnum in steps):
                w.write(i.group(0)[1:-1])
                w.write('\n#Powered by MCPoly\n')
            step = step + 1
    w.close()