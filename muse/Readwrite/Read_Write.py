"""
  MUSE -- A Multi-algorithm-collaborative Universal Structure-prediction Environment

  Copyright (C) 2010-2017 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation
  version 2 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY 
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  Modified form ase.
  E-mail: zl.liu@163.com
"""

import os,math
import numpy as np
from muse.Similarity import Calc_Angle
from muse.Genstrs import Generate_Supercell

def get_atomtypes(fname):
    atomtypes=[]
    if fname.find('.gz') != -1:
        import gzip
        f = gzip.open(fname)
    elif fname.find('.bz2') != -1:
        import bz2
        f = bz2.BZ2File(fname)
    else:
        f = open(fname)
    for line in f:
        if line.find('TITEL') != -1:
            atomtypes.append(line.split()[3].split('_')[0].split('.')[0])
    return atomtypes

def atomtypes_outpot(posfname, numsyms):
    import os.path as op
    import glob

    fnames = [posfname.replace('POSCAR', 'POTCAR').replace('CONTCAR', 
                                                           'POTCAR')]
    fnames.append(posfname.replace('POSCAR', 'OUTCAR').replace('CONTCAR',
                                                               'OUTCAR'))
    # Try the same but with compressed files
    fsc = []
    for fn in fnames:
        fsc.append(fn + '.gz')
        fsc.append(fn + '.bz2')
    for f in fsc:
        fnames.append(f)

    vaspdir = op.dirname(posfname)
    fs = glob.glob(vaspdir + '*POTCAR*')
    for f in fs:
        fnames.append(f)
    fs = glob.glob(vaspdir + '*OUTCAR*')
    for f in fs:
        fnames.append(f)

    tried = []
    files_in_dir = os.listdir('.')
    for fn in fnames:
        if fn in files_in_dir:
            tried.append(fn)
            at = get_atomtypes(fn)
            if len(at) == numsyms:
                return at

    raise IOError('Could not determine chemical symbols. Tried files ' 
                  + str(tried))

def get_atomtypes_from_formula(formula):
    symbols = Generate_Supercell.string2symbols(formula.split('_')[0])
    atomtypes = [symbols[0]]
    for s in symbols[1:]:
        if s != atomtypes[-1]: atomtypes.append(s)
    return atomtypes

def read_vasp(filename='CONTCAR'): 
    if isinstance(filename, str):
        f = open(filename)
    else:
        f = filename

    line1 = f.readline()

    lattice_constant = float(f.readline().split()[0])

    a = []
    for ii in range(3):
        s = f.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        a.append(floatvect)

    basis_vectors = np.array(a) * lattice_constant

    atom_symbols = []
    numofatoms = f.readline().split()

    vasp5 = False
    try:
        int(numofatoms[0])
    except ValueError:
        vasp5 = True
        atomtypes = numofatoms
        numofatoms = f.readline().split()

    # check for comments in numofatoms line and get rid of them if necessary
    commentcheck = np.array(['!' in s for s in numofatoms])
    if commentcheck.any():
        # only keep the elements up to the first including a '!':
        numofatoms = numofatoms[:np.arange(len(numofatoms))[commentcheck][0]]

    if not vasp5:
        atomtypes = line1.split()
       
        numsyms = len(numofatoms)
        if len(atomtypes) < numsyms:
            # First line in POSCAR/CONTCAR didn't contain enough symbols.

            # Sometimes the first line in POSCAR/CONTCAR is of the form
            # "CoP3_In-3.pos". Check for this case and extract atom types
            if len(atomtypes) == 1 and '_' in atomtypes[0]:
                atomtypes = get_atomtypes_from_formula(atomtypes[0])
            else:
                atomtypes = atomtypes_outpot(f.name, numsyms)
        else:
            try:
                for atype in atomtypes[:numsyms]:
                    if not atype in Generate_Supercell.ChemicalSymbols:
                        raise KeyError
            except KeyError:
                atomtypes = atomtypes_outpot(f.name, numsyms)

    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        [atom_symbols.append(atomtypes[i]) for na in xrange(numofatoms[i])]

    # Check if Selective dynamics is switched on
    sdyn = f.readline()
    selective_dynamics = sdyn[0].lower() == "s"

    # Check if atom coordinates are cartesian or direct
    if selective_dynamics:
        ac_type = f.readline()
    else:
        ac_type = sdyn
    cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
    tot_natoms = sum(numofatoms)
    atoms_pos = np.empty((tot_natoms, 3))
    if selective_dynamics:
        selective_flags = np.empty((tot_natoms, 3), dtype=bool)
    for atom in xrange(tot_natoms):
        ac = f.readline().split()
        atoms_pos[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        if selective_dynamics:
            curflag = []
            for flag in ac[3:6]:
                curflag.append(flag == 'F')
            selective_flags[atom] = curflag
    # Done with all reading
    if type(filename) == str:
        f.close()
    if cartesian:
        atoms_pos *= lattice_constant
    atoms = Generate_Supercell.Generate_Supercell(symbols = atom_symbols, cell = basis_vectors, pbc = True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    return atoms

def write_vasp(filename, atoms, label='', direct=False, sort=None, symbol_count = None, long_format=True, vasp5=True):
    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename
    
    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to VASP input")
        else:
            atoms = atoms[0]

    # Write atom positions in scaled or cartesian coordinates
    if direct:
        coord = atoms.get_scaled_positions()
    else:
        coord = atoms.get_positions()

    if sort:
        ind = np.argsort(atoms.get_chemical_symbols())
        symbols = np.array(atoms.get_chemical_symbols())[ind]
        coord = coord[ind]

    else:
        symbols = atoms.get_chemical_symbols()

    # Create a list sc of (symbol, count) pairs
    if symbol_count:
        sc = symbol_count
    else:
        sc = []
        psym = symbols[0]
        count = 0
        for sym in symbols:
            if sym != psym:
                sc.append((psym, count))
                psym = sym
                count = 1
            else:
                count += 1
        sc.append((psym, count))

    # Create the label
    if label == '':
        for sym, c in sc:
            label += '%2s ' % sym
    f.write(label + '\n')

    f.write('%19.16f\n' % 1.0)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write(latt_form % el)
        f.write('\n')

    # If we're writing a VASP 5.x format POSCAR file, write out the
    # atomic symbols
    if vasp5:
        for sym, c in sc:
            f.write(' %3s' % sym)
        f.write('\n')

    # Numbers of each atom
    for sym, count in sc:
        f.write(' %3i' % count)
    f.write('\n')

    #if atoms.constraints:
        #f.write('Selective dynamics\n')

    if direct:
        f.write('Direct\n')
    else:
        f.write('Cartesian\n')

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'
    for iatom, atom in enumerate(coord):
        for dcoord in atom:
            f.write(cform % dcoord)

        f.write('\n')

    if type(filename) == str:
        f.close()


def read_siesta(fname,nu):
    from muse.Readwrite.ReadInput import indict
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
    
    atomorder.sort()
    atomlist = [str(atomorder[i][0]) + str(atomorder[i][1][0]) for i in range(len(atomorder))]
    
    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])

    f = open(fname, 'r')

    cellread = []
    for i in range(3):
        cellread.append([float(x) for x in f.readline().split()])

    natoms = int(f.readline())

    scl_pos = []
    for _ in range(natoms):
        lineread = f.readline().split()[2:]
        scl_pos.append([float(x) for x in lineread])

    atoms = Generate_Supercell.Generate_Supercell(sys_name, cell=cellread, scaled_positions=scl_pos)

    return atoms

def write_siesta(fname,atoms):
    os.system("rm %s"%fname)

    for i in range(3):
        with open(fname,'a') as structfile: print >>structfile,"%21.16f  %21.16f  %21.16f"\
             %(atoms.get_cell()[i][0], atoms.get_cell()[i][1],atoms.get_cell()[i][2])

    with open(fname,'a') as structfile: print >>structfile,len(atoms.get_positions())

    i=1; j=0
    for num in atoms.get_atomic_numbers():
        if j > 0:
            if num != atoms.get_atomic_numbers()[j-1]: i +=1

        with open(fname,'a') as structfile: print >>structfile,"%2d %2d %19.16f  %19.16f  %19.16f"\
             %(i,num,atoms.get_scaled_positions()[j][0], atoms.get_scaled_positions()[j][1],atoms.get_scaled_positions()[j][2])
        j +=1


def read_pwscf(fname,nu):
    from muse.Readwrite.ReadInput import indict
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
    
    atomorder.sort()
    atomlist = [str(atomorder[i][0]) + str(atomorder[i][1][0]) for i in range(len(atomorder))]
    
    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])

    f = open(fname, 'r')

    line = f.readline()

    while line:
        if "CELL_PARAMETERS" in line:
            cellread = []
            scl_pos = []
            for _ in range(3):
                cellread.append([float(x)*0.529177 for x in f.readline().split()])
                
            line = f.readline()
            line = f.readline()
            
            for _ in range(total_atoms):
                line = f.readline()
                lineread = line.split()[1:4]
                #print line
                scl_pos.append([float(x) for x in lineread])
                
        line = f.readline()
        
    #print cellread,scl_pos,total_atoms
    
    atoms = Generate_Supercell.Generate_Supercell(sys_name, cell=cellread, scaled_positions=scl_pos)

    return atoms

def write_pwscf(fname,atoms):
    from muse.Readwrite.ReadInput import indict
    with open(fname,'w') as structfile: print >>structfile,"ATOMIC_SPECIES"
    for x in indict['NameOfAtoms']:
        with open(fname,'a') as structfile: print >>structfile,"%2s  1.0  %2s.UPF"%(x,x)
    
    with open(fname,'a') as structfile: print >>structfile,"CELL_PARAMETERS"
    for i in range(3):
        with open(fname,'a') as structfile: print >>structfile,"%21.16f  %21.16f  %21.16f"\
             %(atoms.get_cell()[i][0]/0.529177, atoms.get_cell()[i][1]/0.529177,atoms.get_cell()[i][2]/0.529177)

    i=1; j=0
    with open(fname,'a') as structfile: print >>structfile,"ATOMIC_POSITIONS {crystal}"
    for num in atoms.get_atomic_numbers():
        if j > 0:
            if num != atoms.get_atomic_numbers()[j-1]: i +=1

        with open(fname,'a') as structfile: print >>structfile,"%2s %19.16f  %19.16f  %19.16f 0 0 0"\
             %(Generate_Supercell.ChemicalSymbols[num],atoms.get_scaled_positions()[j][0], atoms.get_scaled_positions()[j][1],atoms.get_scaled_positions()[j][2])
        j +=1


def read_lammps(fname,nu):
    from muse.Readwrite.ReadInput import indict
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
    
    atomorder.sort()
    atomlist = [str(atomorder[i][0]) + str(atomorder[i][1][0]) for i in range(len(atomorder))]
    
    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])

    f = open(fname, 'r')

    line = f.readline()
    lo = [] ; hi = [] ; tilt = []; scl_pos0 = []
    while line:
        if 'ITEM: BOX BOUNDS' in line:
            
            tilt_items = line.split()[3:]
            for i in range(3):
                line = f.readline()
                fields = line.split()
                lo.append(float(fields[0]))
                hi.append(float(fields[1]))
                if (len(fields) >= 3):
                    tilt.append(float(fields[2]))
                    
            if (len(tilt) >= 3):
                if (len(tilt_items) >= 3):
                    xy = tilt[tilt_items.index('xy')]
                    xz = tilt[tilt_items.index('xz')]
                    yz = tilt[tilt_items.index('yz')]
                    
                else:
                    xy = tilt[0]
                    xz = tilt[1]
                    yz = tilt[2]
                    
            else:
                xy = xz = yz = 0
            xhilo = (hi[0] - lo[0]) - (xy**2)**0.5 - (xz**2)**0.5
            yhilo = (hi[1] - lo[1]) - (yz**2)**0.5
            zhilo = (hi[2] - lo[2])
            if xy < 0:
                if xz < 0:
                    celldispx = lo[0] - xy -xz
                else:
                    celldispx = lo[0] - xy
            else:
                celldispx = lo[0]
            celldispy = lo[1] 
            celldispz = lo[2]
            
            cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]
            celldisp = [[celldispx, celldispy, celldispz]]
            
            line = f.readline()
            for _ in range(total_atoms):
                line = f.readline()
                lineread = line.split()[2:5]
                #print line
                scl_pos0.append([float(x) for x in lineread])
                
        line = f.readline()
    #print scl_pos0

    atoms = Generate_Supercell.Generate_Supercell(sys_name, cell=cell, scaled_positions=scl_pos0)

    return atoms

def write_lammps(fname,atoms):
    from muse.Readwrite.ReadInput import indict
    with open(fname,'w') as structfile: print >>structfile,"Comment line"

    with open(fname,'a') as structfile: print >>structfile,"%d  atoms"%(len(atoms.get_positions()))
    with open(fname,'a') as structfile: print >>structfile,"%d  atom types"%(len(indict['Num_atomsInChemForm']))
    
    cell = atoms.get_cell()
    Len_abc=[]
    for i in range(3):
        Len_abc.append((cell[i][0]**2+cell[i][1]**2+cell[i][2]**2)**.5)
            
    alpha = Calc_Angle.calc_angle(cell[1],cell[2]) * math.pi/180
    beta = Calc_Angle.calc_angle(cell[2],cell[0]) * math.pi/180
    gamma = Calc_Angle.calc_angle(cell[0],cell[1]) * math.pi/180
    
    #print alpha,gamma,beta
    xy= Len_abc[1] * math.cos(gamma)
    xz= Len_abc[2] * math.cos(beta)
    yz= (Len_abc[1]*Len_abc[2]*math.cos(alpha)-xy*xz)/(math.sqrt(Len_abc[1]**2 - xy**2))
    
    newvector = np.matrix([[Len_abc[0],0.0,0.0],\
                           [xy,Len_abc[1],0.0],\
                           [xz,yz,Len_abc[2]]])
    
    with open(fname,'a') as structfile: print >>structfile,"0.000000  %10.6f   xlo xhi"%Len_abc[0]
    with open(fname,'a') as structfile: print >>structfile,"0.000000  %10.6f   ylo yhi"%Len_abc[1]
    with open(fname,'a') as structfile: print >>structfile,"0.000000  %10.6f   zlo zhi"%Len_abc[2]
    with open(fname,'a') as structfile: print >>structfile,"%10.6f  %10.6f  %10.6f   xy xz yz"%(xy,xz,yz)
    with open(fname,'a') as structfile: print >>structfile,"Atoms\n"
    
    i=1; j=0
    for num in atoms.get_atomic_numbers():
        if j > 0:
            if num != atoms.get_atomic_numbers()[j-1]: i +=1
            
        newcoord = atoms.get_scaled_positions()[j]*newvector
        #print list(np.array(newcoord))[0][0]
        with open(fname,'a') as structfile: print >>structfile,"%2d %2d %19.16f  %19.16f  %19.16f"\
                      %(j+1,i,list(np.array(newcoord))[0][0],list(np.array(newcoord))[0][1],list(np.array(newcoord))[0][2])
        j +=1

def read_castep(fname,nu):
    from muse.Readwrite.ReadInput import indict
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
    
    atomorder.sort()
    atomlist = [str(atomorder[i][0]) + str(atomorder[i][1][0]) for i in range(len(atomorder))]
    
    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])
        
    f = open(fname, 'r')
    
    line = f.readline()
    
    while line:
        if "Unit Cell" in line:
            line = f.readline()
            while line:
                if "Unit Cell" in line:
                    cellread = []
                    scl_pos = []
                    chem_symbols = []
                    line = f.readline()
                    line = f.readline()
                    for _ in range(3):
                        cellread.append([float(x) for x in f.readline().split()[:3]])
                    #print cellread
                    
                    for _ in range(16):
                        line = f.readline()
                        #print line
                        
                    for _ in range(total_atoms):
                        line = f.readline()
                        #print line
                        lineread = line.split()[:6]
                        chem_symbols.append(line.split()[1])
                        scl_pos.append([float(x) for x in lineread[3:6]])
                        
                line = f.readline()
        line = f.readline()
        
    atoms = Generate_Supercell.Generate_Supercell(sys_name, cell=cellread, scaled_positions=scl_pos)
    atoms.set_chemical_symbols(chem_symbols)
    
    return atoms

def write_castep(fname,atoms):

    with open(fname,'w') as structfile: print >>structfile,"%BLOCK LATTICE_CART"

    for i in range(3):
        with open(fname,'a') as structfile: print >>structfile,"%21.16f  %21.16f  %21.16f"\
             %(atoms.get_cell()[i][0], atoms.get_cell()[i][1],atoms.get_cell()[i][2])

    with open(fname,'a') as structfile: print >>structfile,"%ENDBLOCK LATTICE_CART"
    
    with open(fname,'a') as structfile: print >>structfile,"%BLOCK POSITIONS_FRAC"

    i=1; j=0
    for num in atoms.get_atomic_numbers():
        if j > 0:
            if num != atoms.get_atomic_numbers()[j-1]: i +=1

        with open(fname,'a') as structfile: print >>structfile,"%2s %19.16f  %19.16f  %19.16f"\
             %(Generate_Supercell.ChemicalSymbols[num],atoms.get_scaled_positions()[j][0], atoms.get_scaled_positions()[j][1],atoms.get_scaled_positions()[j][2])
        j +=1

    with open(fname,'a') as structfile: print >>structfile,"%ENDBLOCK POSITIONS_FRAC"
    
    
def read_abinit(fname,nu):
    from muse.Readwrite.ReadInput import indict
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
    
    atomorder.sort()
    atomlist = [str(atomorder[i][0]) + str(atomorder[i][1][0]) for i in range(len(atomorder))]
    
    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])
        
    f = open(fname, 'r')
    
    line = f.readline()
    
    while line:
        if "(rprimd)" in line.split():
            cellread = []
            for _ in range(3):
                cellread.append([float(x)*.529177 for x in f.readline().split()[:3]])
        line = f.readline()
        
    f.seek(0)
    line = f.readline()
    while line:
        if "(xred)" in line.split():
            scl_pos = []
            for _ in range(total_atoms):
                scl_pos.append([float(x) for x in f.readline().split()[:3]])
        line = f.readline()
        
    atoms = Generate_Supercell.Generate_Supercell(sys_name, cell=cellread, scaled_positions=scl_pos)
    #atoms.set_chemical_symbols(chem_symbols)
    return atoms

def write_abinit(fname,atoms):
    from muse.Readwrite.ReadInput import indict
    with open(fname,'w') as structfile: print >>structfile,"ntypat%6d"%(len(indict['NameOfAtoms']))
    with open(fname,'a') as structfile: print >>structfile,"natom %6d"%(len(atoms.get_positions()))
    
    znucllist = ''
    for x in atoms.get_atomic_numbers():
        znucllist += ("%4s"%str(x))
        
    with open(fname,'a') as structfile: print >>structfile,"znucl    %s"%(znucllist)
    
    i=1; j=0
    typelist = ''
    for num in atoms.get_atomic_numbers():
        if j > 0:
            if num != atoms.get_atomic_numbers()[j-1]:
                i +=1
        typelist += ("%4s"%str(i))
        j +=1
        
    with open(fname,'a') as structfile: print >>structfile,"typat   %s"%(typelist)
    
    with open(fname,'a') as structfile: print >>structfile,"acell      3*1.0"
    
    with open(fname,'a') as structfile: print >>structfile,"rprim   %21.16f  %21.16f  %21.16f"\
         %(atoms.get_cell()[0][0]/0.529177, atoms.get_cell()[0][1]/0.529177,atoms.get_cell()[0][2]/0.529177)
    
    for i in range(2):
        with open(fname,'a') as structfile: print >>structfile,"%29.16f  %21.16f  %21.16f"\
             %(atoms.get_cell()[i+1][0]/0.529177, atoms.get_cell()[i+1][1]/0.529177,atoms.get_cell()[i+1][2]/0.529177)

    with open(fname,'a') as structfile: print >>structfile,"xred    %21.16f  %21.16f  %21.16f"\
         %(atoms.get_scaled_positions()[0][0], atoms.get_scaled_positions()[0][1],atoms.get_scaled_positions()[0][2])
    
    for j in range(len(atoms.get_positions())-1):
        with open(fname,'a') as structfile: print >>structfile,"%29.16f  %21.16f  %21.16f"\
             %(atoms.get_scaled_positions()[j+1][0], atoms.get_scaled_positions()[j+1][1],atoms.get_scaled_positions()[j+1][2])


    
if __name__ == '__main__':
    #atm = read_siesta("mgsio3.STRUCT_IN",4)
    #write_siesta("mgsio3.STRUCT_IN2",atm)
    #print atm
    #atm = read_pwscf("OUT",2)
    
    #write_pwscf("l.in",atm)
    #atm = read_lammps("relax.lammpstrj",26)
    #atm = read_vasp("p2.vasp")
    #write_lammps("data",atm)

    #atm2 = read_lammps("relax.lammpstrj",18)
    #write_vasp("CONTCAR",atm,label='xx', direct=True,sort=True,vasp5=True)
    

    #atm = read_castep("castep.castep",2)
    #write_vasp("CONTCAR",atm,label='xx', direct=True,sort=True,vasp5=True)
    #write_castep("cast.cell",atm)
    atm = read_abinit("abinit.outF",1)
    print atm
    write_vasp("CONTCAR",atm,label='xx', direct=True,sort=True,vasp5=True)
    write_abinit("xcell",atm)
