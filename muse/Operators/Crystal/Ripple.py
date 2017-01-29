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

from numpy import pi,cos
import numpy as np
import os,random,copy,pickle,Mut_Gauss
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict

class RIPPLE:
    def __init__(self,nu,ng,nn,Old_cry,BigDict,V_best,picked_indi,Num_OperSuccess,kT,dr=0.3):
        self.nu = nu
        self.ng = ng
        self.nn = nn
        self.Old_cry = Old_cry
        self.V_best = V_best
        self.BigDict = BigDict
        self.picked_indi = picked_indi
        self.Num_OperSuccess = Num_OperSuccess
        self.kT = kT
        self.dr = dr

    def Ripple(self,O_cry):
        Old_cry = copy.deepcopy(O_cry)
        pos = Old_cry.get_positions()
        ax = random.sample(range(3),3)

        for i in range(len(pos)):
            pos[i][ax[0]] += self.dr*cos(2*pi*1*pos[i][ax[1]] + 2*pi*random.random())*\
                                       cos(2*pi*1*pos[i][ax[2]] + 2*pi*random.random())

        Old_cry.set_positions(pos)

        return Old_cry

    def Run_Ripple(self):
        O_cry = self.Ripple(self.Old_cry[self.picked_indi][1])
        spg0 = Findspg.Findspg(O_cry)
        with open('../log.muse','a') as logfile:
            print >>logfile, "Mutate and ripple for %02d: %s %s"%(self.picked_indi+1,spg0[0],spg0[1])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            wlabel = indict['NameSys'][0]+": "+str(self.ng)+'-'+str(self.picked_indi+1)
            Read_Write.write_vasp('POSCAR',O_cry,label=wlabel,direct=True,sort=True,vasp5=True)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile: print >>logfile, self.Old_cry[self.picked_indi][0],enth,spg1[0]

            if enth <= self.Old_cry[self.picked_indi][0]:
                self.Num_OperSuccess[self.ng][5] +=1
                ifdo = 0

            else:
                ifdo -= 1
                if indict['Algorithm'][0] == 'MAC' or indict['Algorithm'][0] == 'SA':
                    if np.exp((self.Old_cry[self.picked_indi][0] - enth) / self.kT) < np.random.uniform():
                        BigDict[self.nu][self.ng][self.Old_cry[self.picked_indi][0]] = self.Old_cry[self.picked_indi][1]

                    #if enth in BigDict[nu][ng].keys():
                        #del BigDict[nu][ng][enth]

            if 111111.11111111 in BigDict[self.nu][self.ng].keys(): del self.BigDict[self.nu][self.ng][111111.11111111]

        with open('.Num_OperSuccess0','w') as Num_OperSuccessFile:
            pickle.dump(self.Num_OperSuccess,Num_OperSuccessFile)

        os.system("mv .Num_OperSuccess0 .Num_OperSuccess")

        return nk,self.Num_OperSuccess,BigDict

    def Run_MutRipple(self):
        #GaussStdDev = float(indict['InitGaussStdDev'][0])
        mutgauss = Mut_Gauss.MUTGAUSS(self.nu,self.ng,self.nn,self.Old_cry[self.picked_indi][1],self.BigDict,self.picked_indi,self.V_best,self.Num_OperSuccess,self.kT)
        O_cry = mutgauss.Mut_Gauss()
        O_cry = self.Ripple(O_cry)
        spg0 = Findspg.Findspg(O_cry)
        with open('../log.muse','a') as logfile: print >>logfile, "Ripple for %02d: %s"%(self.picked_indi+1,spg0[0])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            wlabel = indict['NameSys'][0]+": "+str(self.ng)+'-'+str(self.picked_indi+1)
            Read_Write.write_vasp('POSCAR',O_cry,label=wlabel,direct=True,sort=True,vasp5=True)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile: print >>logfile, self.Old_cry[self.picked_indi][0],enth,spg1[0]

            if enth <= self.Old_cry[self.picked_indi][0]:
                self.Num_OperSuccess[self.ng][6] +=1
                ifdo = 0

            else:
                ifdo -= 1
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
