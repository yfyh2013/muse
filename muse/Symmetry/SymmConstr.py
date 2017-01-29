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

import os
import numpy as np
from muse.Symmetry import Findspg
from muse.Crystaloper import Splitnum,Balance,HardConstrain
from muse.Readwrite import Read_Write
from muse.Genstrs import Generate_Supercell,Run_genpc
from spglib import spglib
from muse.Readwrite.ReadInput import indict


def GenCryAndReduce(nu,sys_name,Num_atomsInChemForm,sg_dict):
    spg = 0
    maxtry = 0
    while spg == 1 or spg in sg_dict:
        maxtry += 1
        #print "maxtry = %d"%maxtry
        randstr = []

        ntry = 0
        while True:
            ntry +=1
            randstr0,sg_num,Isdirectmatch = Run_genpc.Run_genpc(nu,sg_dict)

            randstr = randstr0
            if Isdirectmatch: break

            if not Isdirectmatch:
                nnn = Splitnum.Splitnum(nu)

                n1,n2,n3 = Balance.Balance(randstr0,nnn)
                #print nnn,n1,n2,n3

                if nu*sum(Num_atomsInChemForm) == len(randstr0.get_positions()):
                    randstr = randstr0
                    break
                elif sum(Num_atomsInChemForm) == len(randstr0.get_positions()):
                    randstr = randstr0.repeat((n1,n2,n3))
                    break
                elif nu*sum(Num_atomsInChemForm) % len(randstr0.get_positions()) == 0:
                    mmm = Splitnum.Splitnum(int(nu*sum(Num_atomsInChemForm) / len(randstr0.get_positions())))
                    m1,m2,m3 = Balance.Balance(randstr0,mmm)
                    #print mmm,m1,m2,m3
                    randstr = randstr0.repeat((m1,m2,m3))
                    break

                pc = spglib.find_primitive(randstr0,symprec=0.01,angle_tolerance=0.01)

                if pc[0] is not None:
                    if sum(Num_atomsInChemForm) == len(pc[1]):
                        pcell = []
                        for i in range(len(pc[0])):
                            pcell.append([pc[0][i][j] for j in range(3)])

                        V_init = abs(np.linalg.det(pcell))
                        for ll in range(len(pcell)):
                            for jj in range(3):
                                pcell[ll][jj] *= (nu*float(indict['V_guess'][0])/V_init)**(1./3.)

                        pos = []
                        for i in range(len(pc[1])):
                            pos.append([pc[1][i][j] for j in range(3)])

                        randstr0= Generate_Supercell.Generate_Supercell(sys_name,cell=pcell,scaled_positions=pos)
                        randstr = randstr0.repeat((n1,n2,n3))

                        break


            if ntry > 500: break

        randstr.set_cell(randstr.get_cell()*(nu*float(indict['V_guess'][0])/randstr.get_volume())**(1./3.),scale_atoms=True)

        #print len(randstr.get_scaled_positions())
        #print len(set(list(randstr.get_scaled_positions())))
        # if randstr:
        #    prec = float(indict['SymPrec'][0])
        #    spg0 = spglib.get_spacegroup(randstr, symprec=prec).split()
        #    spg = spg0[1][1:-1]
        #print spg
        #spg = sg_num

        if maxtry >= 100:
            #print "%d used anyway!"%spg
            break

    return randstr,sg_num

def SymmConstr(nu,ng,nn,sg_dict):
    with open('../log.muse','a') as logfile: print >>logfile, "Random str %02d..."%nn

    ifempty = False
    Num_atomsInChemForm = [int(i) for i in indict['Num_atomsInChemForm']]

    atomorder = {indict['NameOfAtoms'][i]:\
                [int(indict['Num_atomsInChemForm'][i]),\
                 int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()

    atomorder.sort()

    atomlist = [str(atomorder[ii][0]) + str(atomorder[ii][1][0]) for ii in range(len(atomorder))]

    sys_name = ''
    total_atoms=0
    for i in range(len(atomlist)):
        sys_name = sys_name + atomlist[i]
        total_atoms = total_atoms + int(atomorder[i][1][0])

    randstr,sg_num = GenCryAndReduce(nu,sys_name,Num_atomsInChemForm,sg_dict)

#
#    L_len=True; L_dist=True; L_angle=True
#    mtry = 1
#    while L_len or L_dist or L_angle:
#        randstr,sg_num = GenCryAndReduce(nu,sys_name,Num_atomsInChemForm,sg_dict)
#
#        if len(randstr) != 0:
#            L_dangle,L_angle,L_len,L_dist = HardConstrain.HardConstrain(randstr,nu)
#        if mtry > 100:
#            ifempty = True
#            break
#        mtry +=1


    prec = float(indict['SymPrec'][0])
    spg00 = spglib.get_spacegroup(randstr, symprec=prec)
    if spg00:
        spg0 = spg00.split()
        spg1 = [spg0[0], spg0[1][1:-1]]

        wlabel = indict['NameSys'][0]+": "+"%02d-%02d"%(ng,nn)+'  '+"Random Str.: "+spg1[0]+spg1[1]
        Read_Write.write_vasp('POSCAR',randstr,label=wlabel,direct=True,sort=True,vasp5=True)
        os.system("cat POSCAR >> ../poscars_notopt")


        #print "Created crystal is saved: %s"%spg1[1]
        sg_dict.append(int(spg1[1]))
        #sg_dict.append(sg_num)
        #sg_dict.append(spg1[1]])

    return sg_dict,ifempty

if __name__ == '__main__':
    sg_dict = []
    for nn in range(30):
        sg_dict,ifempty = SymmConstr(4,1,nn,sg_dict)
        if ifempty: sg_dict = []; print 'sg_dict emptied.'
        print sg_dict

    sg_dict.sort()
    print sg_dict
