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
import os,random,copy,pickle
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Crystaloper import HardConstrain
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict

class MUTGAUSS:
    def __init__(self,nu,ng,nn,Old_cry,BigDict,picked_indi,V_best,Num_OperSuccess,kT):
        self.nu = nu
        self.ng = ng
        self.nn = nn
        self.BigDict = BigDict
        self.Old_cry = Old_cry
        self.V_best = V_best
        self.picked_indi = picked_indi
        self.GaussStdDev = float(indict['InitGaussStdDev'][0])
        self.Num_OperSuccess = Num_OperSuccess
        self.kT = kT

    def Mut_Gauss(self,supercell):
        Max_try = 0
        L_len=True; L_dist=True; L_angle=True; L_dangle=True
        while L_len or L_dist or L_angle:
            #print self.Old_cry
            O_cry = copy.deepcopy(supercell)
            ee=[]
            for i in range(6):
                ee.append(random.gauss(0.0,self.GaussStdDev))

            mut_matrix=np.matrix([[1.0+ee[0],0.5*ee[1],0.5*ee[2]],\
                                  [0.5*ee[1],1.0+ee[3],0.5*ee[4]],\
                                  [0.5*ee[2],0.5*ee[4],1.0+ee[5]]])

            O_cry.set_cell(mut_matrix*O_cry.get_cell(),scale_atoms=True)

            O_cry.set_cell(O_cry.get_cell()*(self.V_best/O_cry.get_volume())**(1./3.),scale_atoms=True)

            L_dangle,L_angle,L_len,L_dist = HardConstrain.HardConstrain(O_cry,self.nu)

            Max_try +=1
            if Max_try >= 1000: break

        Read_Write.write_vasp('POSCAR',O_cry,label=indict['NameSys'][0]+": "+str(self.ng)+'-'+str(self.nn),direct=True,sort=True,vasp5=True)

        return O_cry

    def Run_Mut(self):
        spg0 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Mutate %.3f for %02d: %s"%(self.GaussStdDev,self.picked_indi+1,spg0[0])

            Old_cry = self.Mut_Gauss(self.Old_cry[self.picked_indi][1])
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,Old_cry)
            with open('../log.muse','a') as logfile: print >>logfile, self.Old_cry[self.picked_indi][0],enth

            #spg1 = Findspg.Findspg(BigDict[nu][ng][enth])

            if enth <= self.Old_cry[self.picked_indi][0]:# or spg1[0] != spg0[0]:
                self.Num_OperSuccess[self.ng][0] +=1
                ifdo = 0

            else:
                ifdo -= 1
                self.GaussStdDev -= 0.1
                if indict['Algorithm'][0] == 'MAC' or indict['Algorithm'][0] == 'SA':
                    if np.exp((self.Old_cry[self.picked_indi][0] - enth) / self.kT) < np.random.uniform():
                        BigDict[self.nu][self.ng][self.Old_cry[self.picked_indi][0]] = self.Old_cry[self.picked_indi][1]

                    #if enth in BigDict[nu][ng].keys():
                        #del BigDict[nu][ng][enth]

            if 111111.11111111 in BigDict[self.nu][self.ng].keys(): del BigDict[self.nu][self.ng][111111.11111111]

        with open('.Num_OperSuccess0','w') as Num_OperSuccessFile:
            pickle.dump(self.Num_OperSuccess,Num_OperSuccessFile)

        os.system("mv .Num_OperSuccess0 .Num_OperSuccess")

        return nk,self.Num_OperSuccess,BigDict
