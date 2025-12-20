import os
import re
import time as t
import shutil

def pltvib(file, nums, orcaloc='./', fileloc='./'):
    """
    The method to get the movement of some frequencies about .hess file by ORCA.
    pltvib(file, nums=[], orcaloc='./', fileloc='./')
    file: Your .inp file name.
    nums: The frequency numbers.
    orcaloc: Your location of ORCA.
             If you have set ORCA from https://www.orcasoftware.de/tutorials_orca/first_steps/trouble_install.html#path-variable, 
             you can skip this key word.
    fileloc: Your location of .hess file. The default is your current location.
    saveloc: Your location of .inp file. The default is your current location.
    Example:
        Input:
            from MCPoly.orcaset import XYZtoINP
            from MCPoly.orcaset import orca
            
            XYZtoINP('Atom2', method='B3LYP', basis_set='def2-TZVP', opt=True, freq=True, \
                        maxcore=4096, corenum=8)
            orca('Atom2', orcaloc='./MCPoly/orca/')
        Output:
            <<COMMAND LINES>>
            Program Completed.[Mon Apr 17 14:25:12 2023] -> See the current location to check files.
            
    TIPS: 1. Don't put the .inp to your ORCA installment file,  it may disturb your work. 
          2. If you don't input key word 'orcaloc' while you don't set ORCA,  you will only get an empty .out file in the document.
             There will be a warning: sh: orca: 'command not found'.
          3. If you have set key word 'maxcore' and 'corenum',  you must input 'orcaloc'.
             See in https://www.orcasoftware.de/tutorials_orca/first_steps/parallel.html
          4. You can see the error of ORCA directly from the output.
    """

    for num in nums:
        if orcaloc == './':
            os.system('orca_pltvib {0}{1}.hess {2}'.format(fileloc, file, num))
        else:
            path = os.getcwd() + '/'
            os.chdir(orcaloc)
            print(os.getcwd())
            if fileloc == './':
                fileloc = path
            os.system('./orca_pltvib {0}{1}.hess {2}'.format(fileloc, file, num))
            print('Program Completed.[{0}]\n'.format(t.ctime(t.time())))
            os.chdir(path[:-1])