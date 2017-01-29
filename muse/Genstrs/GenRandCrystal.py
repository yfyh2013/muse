"""
  MUSE -- A Multi-algorithm-collaborative Universal Structure-prediction Environment

  Copyright (C) 2010-2017 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 2 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""

import random,os
import numpy as np
from muse.Readwrite import Read_Write
from muse.Readwrite.ReadInput import indict
from muse.Genstrs import Generate_Supercell
from muse.Crystaloper import HardConstrain,Centering


def GenRandCrystal(nu,ng,nn):
    with open('../log.muse','a') as logfile: print >>logfile, "Random str %02d..."%nn

    global Nopt
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()

    atomorder.sort()
    atomlist = [str(atomorder[ii][0]) + str(atomorder[ii][1][0]) for ii in range(len(atomorder))]

    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])

    L_len=True; L_dist=True; L_angle=True; L_dangle=True

    while L_len or L_dist or L_angle:

        atoms_scl_pos = np.random.uniform(0.0, 1.0, (total_atoms, 3))
        cellp=[]
        for i in range(3): cellp.append(random.random())
        for i in range(3): cellp.append(45 + 90 * random.random())

        V_init = abs(np.linalg.det(Generate_Supercell.cellpar_to_cell(cellp)))
        randstr= Generate_Supercell.Generate_Supercell(sys_name, cell=Generate_Supercell.cellpar_to_cell(cellp)*\
                                                        (nu*float(indict['V_guess'][0])/V_init)**(1./3.),\
                                                        scaled_positions=atoms_scl_pos)

        L_dangle,L_angle,L_len,L_dist = HardConstrain.HardConstrain(randstr,nu)

    randstr = Centering.Centering(randstr)

    wlabel = indict['NameSys'][0]+": "+str(ng)+'-'+str(nn)

    if int(indict['OptCode'][0]) == 1:
        Read_Write.write_vasp('POSCAR',randstr,label=wlabel,direct=True,sort=True,vasp5=True)

    elif int(indict['OptCode'][0]) == 2:
        Read_Write.write_vasp('POSCAR',randstr,label=wlabel,direct=True,sort=True,vasp5=True)
        randstr2 = Read_Write.read_vasp('POSCAR')
        Read_Write.write_siesta("%s.STRUCT_IN"%indict['NameSys'][0],randstr2)
        os.system("rm POSCAR")

    return None
