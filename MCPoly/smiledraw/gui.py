# import numpy as np
import os
import re
import py3Dmol
from ase.io import read, write
from ase import Atoms
from ase.visualize import view
import ipywidgets as widgets
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import AllChem
from ipywidgets import interact
import sys
from PIL import Image

mydir = os.path.dirname( __file__ )
molesdir = os.path.join(mydir, '..', 'moldraw')
sys.path.append(molesdir)
from smiledraw import Smolecule
from smiledraw import Ssub_selection
#from smiledraw import Sbind_selection
from smiledraw import Sring_selection
from smiledraw import S_C_selection

mydir = os.path.dirname( __file__ )
molesdir = os.path.join(mydir, '..', 'moldraw')
sys.path.append(molesdir)
from moldraw import molecule

mydir = os.path.dirname( __file__ )
orcasdir = os.path.join(mydir, '..', 'orcaset')
sys.path.append(orcasdir)
from orcaset import view3dchoose

mydir = os.path.dirname( __file__ )
viewdir = os.path.join(mydir, '..', 'view3d')
sys.path.append(viewdir)
from view3d import view3d

def smi2conf(smiles):
    '''Convert SMILES to rdkit.Mol with 3D coordinates'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        return mol
    else:
        return None

def MolTo3DView(mol, size=(400, 500), surface=False, opacity=0.5):
    """Draw molecule in 3D (Not Origin)
    
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'cartoon'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    mblock = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({'sphere':{'radius':0.4}, 'stick':{'radius':0.1}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    Chem.MolToXYZFile(mol, 'Mole_Untitled.xyz')
    return viewer

def cleanup_qm9_xyz(fname):
    ind = open(fname).readlines()
    nAts = int(ind[0])
    # There are two smiles in the data: the one from GDB and the one assigned from the
    # 3D coordinates in the QM9 paper using OpenBabel (I think).
    gdb_smi, relax_smi = ind[-2].split()[:2]
    ind[1] = '\n'
    ind = ind[:nAts+2]
    for i in range(2, nAts+2):
        l = ind[i]
        l = l.split('\t')
        l.pop(-1)
        ind[i] = '\t'.join(l)+'\n'
    ind = ''.join(ind)
    return ind, gdb_smi, relax_smi

def draw_with_spheres_choose(mol, choose, width, height):
    v = py3Dmol.view(width=width, height=height)
    IPythonConsole.addMolToView(mol, v)
    v.zoomTo()
    v.setStyle({'sphere':{'radius':0.3}, 'stick':{'radius':0.1}})
    for i in choose:
        v.setStyle({'model': -1, 'serial': i},\
                   {'sphere':{'radius':0.5, 'color':'Blue'}, 'stick':{'radius':0.1}})
    v.show()

def view3dchoose(filename, loc='./', choose=[], width=400, height=550):
    path = os.getcwd()
    os.chdir(loc)
    ind = open('{0}.xyz'.format(filename), 'r+')
    mains = ind.readlines()
    main = ''
    for i, line in enumerate(mains):
        if i == 1:
            main = main + '\n'
        else:
            main = main + mains[i]
    #print(main)
    raw_mol = Chem.MolFromXYZBlock(main)
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(conn_mol)
    v = draw_with_spheres_choose(conn_mol, choose, width, height)
    os.chdir(path)
    return v

def Ssub_symbol_mean(substrate):
    #substrates = ['H', 'F', 'Cl', 'Br', 'I', 'OH', 'OMe', 'SH', 'SMe', 'NH2', 'NMe2',\
    #              'Me', 'Et', 'n-Pr', 'i-Pr', 'n-Bu', 't-Bu', 'CF3', 'CCl3', 'NO2', \
    #              'C2H3', 'C=NH', 'CHO', 'COOH', 'COMe', 'CN', 'Ph'\
    #              'cyclo3', 'cyclo5', 'cyclo6', 'SOMe', 'SO2Me', 'SO3H']
    
    if substrate == 'NMe2':
        return 'N(C)C'
    elif substrate == 'Et':
        return 'CC'
    elif substrate == 'n-Pr':
        return 'CCC'
    elif substrate == 'i-Pr':
        return 'C(C)C'
    elif substrate == 'n-Bu':
        return 'CCCC'
    elif substrate == 'i-Bu':
        return 'C(C)CC'
    elif substrate == 't-Bu':
        return 'C(C)(C)C'
    elif substrate == 'C2H3':
        return 'C=C'
    elif substrate == 'COMe':
        return 'C(=O)C'
    elif substrate == 'H':
        return 'H'
    elif 'H' in substrate and substrate not in ('CHO', 'COOH', 'SO3H'):
        a = re.split(r'H', substrate)
        total = ''
        for i in a:
            total = total + i
        return total
    elif 'Me' in substrate[-2:]:
        return substrate[:-2] + 'C'
    else:
        return substrate

def gui():
    """
        The method to draw molecule and save it in the folder.
        You will know how to do as long as you get into the GUI platform.
    """
    
    def show_geos(smiles):
        smiles_detail = Smolecule(smiles)
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def sub_pattern_judge(smiles, number1):
        global command
        pattern_avail = []
        for option in ['1', 'R', 'S', 'Z', 'E', '2']:
            try:
                x = Ssub_selection({option:command}, smiles, number1)
                pattern_avail.append(option)
            except:
                pass
        return pattern_avail

    def begin(button):
        global command
        command = 'H'
        smiles_detail = Smolecule(smiles.value)
        frags_len = len(smiles_detail.frags)
        number1.max = frags_len - 1
        number2.max = frags_len - 1
        main_output1.clear_output()
        with main_output1:
            display(tab_nest)

    def chooseSsub(Ssubbutton):
        global command
        with main_output1:
            pattern_avail = sub_pattern_judge(smiles.value, number1.value)
            pattern.options = pattern_avail
            command = Ssubbutton.description
            command = Ssub_symbol_mean(command)
        #print(command)

    def chooseSring(Sringbutton):
        global command
        with main_output1:
            pattern_avail = sub_pattern_judge(smiles.value, number1.value)
            pattern2.options = pattern_avail
            pattern_avail2 = sub_pattern_judge(smiles.value, number2.value)
            pattern2_2.options = pattern_avail2
        #print(command)

    def chooseSC(SCbutton):
        global command
        with main_output1:
            pattern_avail = sub_pattern_judge(smiles.value, number1.value)
            pattern.options = pattern_avail
            command = SCbutton.description
        #print(command)

    def SsubPremodify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles_detail.S_sub({pattern.value: command}, number1.value)
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def SsubModify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles.value = smiles_detail.S_sub({pattern.value: command}, number1.value)
        frags_len = len(smiles_detail.frags)
        number1.max = frags_len - 1
        number2.max = frags_len - 1
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def SringPremodify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles_detail.S_ring({number1.value: pattern2.value,\
                              number2.value: pattern2_2.value})
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def SringModify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles.value = smiles_detail.S_ring({number1.value: pattern2.value,\
                                             number2.value: pattern2_2.value})
        frags_len = len(smiles_detail.frags)
        number1.max = frags_len - 1
        number2.max = frags_len - 1
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass
    
    def SCPremodify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles_detail.S_C(command, number1.value)
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def SCModify(button):
        global command
        smiles_detail = Smolecule(smiles.value)
        smiles.value = smiles_detail.S_C(command, number1.value)
        frags_len = len(smiles_detail.frags)
        number1.max = frags_len - 1
        number2.max = frags_len - 1
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            print(smiles_detail.smiles)
            d = rdMolDraw2D.MolDraw2DCairo(480, 450)
            rdMolDraw2D.PrepareAndDrawMolecule(d, smiles_detail.image)
            d.WriteDrawingText('untitled.png')
            display(Image.open('untitled.png'))
        with main_output3:
            try:
                smiles_detail.view3d(size=(480, 450))
            except:
                pass

    def normal_highlight1(num1):
        pattern_avail = sub_pattern_judge(smiles.value, num1)
        pattern.options = pattern_avail
        mol = Chem.MolFromSmiles(smiles.value)
        patt = mol
        hit_ats = list(mol.GetSubstructMatch(patt))
        hit_bonds = []
        for bond in patt.GetBonds():
           aid1 = hit_ats[bond.GetBeginAtomIdx()]
           aid2 = hit_ats[bond.GetEndAtomIdx()]
           hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        d = rdMolDraw2D.MolDraw2DCairo(480, 450)
        colours = [(1,0.2,0.2),(1,1,1)]
        atom_cols = {}
        for i, at in enumerate(hit_ats):
            if i == num1:
                atom_cols[at] = colours[0]
            else:
                atom_cols[at] = colours[1]
        bond_cols = {}
        for i, bd in enumerate(hit_bonds):
            bond_cols[bd] = colours[1]

        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,\
                                           highlightAtomColors=atom_cols,\
                                           highlightBonds=hit_bonds,\
                                           highlightBondColors=bond_cols)
        d.WriteDrawingText('untitled.png')
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            display(Image.open('untitled.png'))
        with main_output3:
            view3dchoose('Mole_Untitled', loc='./', choose=[num1], width=480, height=450)

    def normal_highlight2(num1, num2):
        pattern_avail = sub_pattern_judge(smiles.value, num1)
        try:
            pattern_avail.remove("2")
        except:
            pass
        pattern2.options = pattern_avail
        pattern_avail2 = sub_pattern_judge(smiles.value, num2)
        try:
            pattern_avail2.remove("2")
        except:
            pass
        pattern2_2.options = pattern_avail2
        mol = Chem.MolFromSmiles(smiles.value)
        patt = mol
        hit_ats = list(mol.GetSubstructMatch(patt))
        hit_bonds = []
        for bond in patt.GetBonds():
           aid1 = hit_ats[bond.GetBeginAtomIdx()]
           aid2 = hit_ats[bond.GetEndAtomIdx()]
           hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        d = rdMolDraw2D.MolDraw2DCairo(480, 450)
        colours = [(1,0.2,0.2),(1,1,1)]
        atom_cols = {}
        for i, at in enumerate(hit_ats):
            if i in (num1, num2):
                atom_cols[at] = colours[0]
            else:
                atom_cols[at] = colours[1]
        bond_cols = {}
        for i, bd in enumerate(hit_bonds):
            bond_cols[bd] = colours[1]

        rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,\
                                           highlightAtomColors=atom_cols,\
                                           highlightBonds=hit_bonds,\
                                           highlightBondColors=bond_cols)
        d.WriteDrawingText('untitled.png')
        main_output2.clear_output()
        main_output3.clear_output()
        with main_output2:
            display(Image.open('untitled.png'))
        with main_output3:
            view3dchoose('Mole_Untitled', loc='./', choose=[num1, num2], width=480, height=450)
            
    #Begin
    main_output1 = widgets.Output()
    main_output2 = widgets.Output()
    main_output3 = widgets.Output()
    smiles = widgets.Text(description='SMART:',\
                          placeholder = 'Input SMART structure')
    
    beginmodify_button = widgets.Button(description='Sure')
    beginmodify_button.on_click(begin)
    with main_output1:
        display(widgets.HBox([smiles, beginmodify_button]))

    S_sub_mode = widgets.Output()
    S_bind_mode = widgets.Output()
    S_ring_mode = widgets.Output()
    S_C_mode = widgets.Output()
    
    number1 = widgets.BoundedIntText(value=0, min=0, max=100, \
                                     description='Number 1:')
    number2 = widgets.BoundedIntText(value=0, min=0, max=100, \
                                     description='Number 2:')
    pattern = widgets.Dropdown(options=['1', 'R', 'S', '2'], value='1',\
                               description='Structure')
    pattern2 = widgets.Dropdown(options=['1', 'R', 'S'], value='1',\
                                description='Structure')
    pattern2_2 = widgets.Dropdown(options=['1', 'R', 'S'], value='1',\
                                  description='Structure')
    site = widgets.BoundedIntText(value=0, min=0, max=5, \
                                     description='Site:')
    
    tab_nest = widgets.Tab()
    tab_nest.children = [S_sub_mode, S_bind_mode, S_ring_mode, S_C_mode]
    tab_nest.titles = ('S_sub_mode', 'S_bind_mode', 'S_ring_mode', 'S_C_mode')
    
    substrates = ['H', 'F', 'Cl', 'Br', 'I', 'OH', 'OMe', 'SH', 'SMe', 'NH2', 'NMe2',\
                  'Me', 'Et', 'n-Pr', 'i-Pr', 'n-Bu', 't-Bu', 'CF3', 'CCl3', 'NO2', \
                  'C2H3', 'C=NH', 'CHO', 'COOH', 'COMe', 'CN', 'Ph', 'NA',\
                  'cyclo3', 'cyclo5', 'cyclo6', 'SOMe', 'SO2Me', 'SO3H']
    Ssubbuttons1 = []
    sub_premodify_button = widgets.Button(description='Browse')
    sub_modify_button = widgets.Button(description='Sure')
    for x in substrates:
        Ssubbutton = widgets.Button(layout = widgets.Layout(height='30px',\
                                                      min_width='120px'), description=x)
        Ssubbuttons1.append(Ssubbutton)
        Ssubbutton.on_click(chooseSsub)
    SsubBox = widgets.HBox(
        children = [Ssubbuttons1[x] for x in range(len(substrates))])
    with S_sub_mode:
        normal_highlight1_out = widgets.interactive_output(normal_highlight1, {'num1': number1})
        display(widgets.VBox([SsubBox, widgets.HBox([number1, pattern2,\
                                                     sub_premodify_button, sub_modify_button])]))
    sub_premodify_button.on_click(SsubPremodify)
    sub_modify_button.on_click(SsubModify)

    ring_premodify_button = widgets.Button(description='Browse')
    ring_modify_button = widgets.Button(description='Sure')
    with S_ring_mode:
        normal_highlight2_out = widgets.interactive_output(normal_highlight2, {'num1': number1,\
                                                                               'num2': number2})
        display(widgets.VBox([widgets.HBox([number1, pattern2]),\
                              widgets.HBox([number2, pattern2_2]),\
                              widgets.HBox([ring_premodify_button, ring_modify_button])]))
    ring_premodify_button.on_click(SringPremodify)
    ring_modify_button.on_click(SringModify)
    
    Cs = ['O', 'N', 'S', 'C(=O)', 'C(=N)', 'S(=O)', 'S(=O)(=O)']
    SCbuttons1 = []
    C_premodify_button = widgets.Button(description='Browse')
    C_modify_button = widgets.Button(description='Sure')
    for x in Cs:
        SCbutton = widgets.Button(layout = widgets.Layout(height='30px',\
                                                      min_width='120px'), description=x)
        SCbuttons1.append(SCbutton)
        SCbutton.on_click(chooseSC)
    SCBox = widgets.HBox(children = [SCbuttons1[x] for x in range(len(Cs))])
    with S_C_mode:
        normal_highlight1_out = widgets.interactive_output(normal_highlight1, {'num1': number1})
        display(widgets.VBox([SCBox, widgets.HBox([number1,\
                                                   C_premodify_button, C_modify_button])]))
    C_premodify_button.on_click(SCPremodify)
    C_modify_button.on_click(SCModify)
    
    
    viewout_begin = widgets.interactive_output(show_geos, {'smiles': smiles})
    
    box_layout1 = widgets.Layout(overflow='scroll hidden',\
                                border='2px solid black', width='97.5%', height='200px', display='flex')
    box_layout2 = widgets.Layout(overflow='scroll hidden',\
                                border='2px solid black', width='48.5%', height='450px', display='flex')
    main_output1.layout = box_layout1
    main_output2.layout = box_layout2
    main_output3.layout = box_layout2
    
    display(widgets.VBox([main_output1, widgets.HBox([main_output2, main_output3])]))