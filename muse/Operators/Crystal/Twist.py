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

import numpy as np
import os,copy,math,random,pickle
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict
from muse.Crystaloper import Centering

class TWIST:
    def __init__(self,nu,ng,nn,O_cry,BigDict,picked_indi,Num_OperSuccess,kT):
        self.nu = nu
        self.ng = ng
        self.nn = nn
        self.O_cry = O_cry
        self.BigDict = BigDict
        self.picked_indi = picked_indi
        self.Num_OperSuccess = Num_OperSuccess
        self.kT = kT
        self.theta = random.random() * math.pi/2

    def Twist(self,atm):
        atoms = copy.deepcopy(atm)
        atoms = Centering.Centering(atoms)
        pos = atoms.get_scaled_positions()

        for i in range(len(pos)):
            rotate_matrix=np.matrix([[math.cos(self.theta),math.sin(self.theta),0.0],\
                                     [-math.sin(self.theta),math.cos(self.theta),0.0],\
                                     [0.0,0.0,1.0]])
            pos[i] = np.matrix(pos[i]) * rotate_matrix
            self.theta += math.pi/8

        atoms.set_scaled_positions(pos)

        return atoms

    def Run_Twist(self):
        Old_cry = copy.deepcopy(self.O_cry)
        spg0 = Findspg.Findspg(Old_cry[self.picked_indi][1])
        with open('../log.muse','a') as logfile: print >>logfile, "Twist for %02d: %s %s"%(self.picked_indi+1,spg0[0],spg0[1])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            slipatoms = self.Twist(Old_cry[self.picked_indi][1])
            wlabel = indict['NameSys'][0]+": "+str(self.ng)+'-'+str(self.picked_indi+1)
            Read_Write.write_vasp('POSCAR',slipatoms,label=wlabel,direct=True,sort=True,vasp5=True)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile: print >>logfile, Old_cry[self.picked_indi][0],enth,spg1[0]

            if enth < Old_cry[self.picked_indi][0]:# or spg1[0] != spg0[0]:
                ifdo = 0
                if indict['Algorithm'][0] == 'EA2':
                    self.Num_OperSuccess[self.ng][3] +=1
                else:
                    self.Num_OperSuccess[self.ng][6] +=1

            else:
                ifdo -= 1
                if indict['Algorithm'][0] == 'MAC' or indict['Algorithm'][0] == 'SA':
                    if np.exp((Old_cry[self.picked_indi][0] - enth) / self.kT) > np.random.uniform():
                        BigDict[self.nu][self.ng][Old_cry[self.picked_indi][0]] = Old_cry[self.picked_indi][1]

                    #if enth in BigDict[nu][ng].keys():
                        #del BigDict[nu][ng][enth]

            if 111111.11111111 in BigDict[self.nu][self.ng].keys(): del BigDict[self.nu][self.ng][111111.11111111]

        with open('.Num_OperSuccess0','w') as Num_OperSuccessFile:
            pickle.dump(self.Num_OperSuccess,Num_OperSuccessFile)

        os.system("mv .Num_OperSuccess0 .Num_OperSuccess")

        return nk,self.Num_OperSuccess,BigDict
