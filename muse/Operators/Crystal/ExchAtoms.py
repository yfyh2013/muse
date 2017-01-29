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

import os,random,copy,pickle,Mut_Gauss
import numpy as np
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict

class EXCHATOMS:
    def __init__(self,nu,ng,nn,Old_cry,BigDict,V_best,picked_indi,picked_indj,Num_OperSuccess,kT):
        self.nu = nu
        self.ng = ng
        self.nn = nn
        self.Old_cry = Old_cry
        self.BigDict = BigDict
        self.V_best = V_best
        self.picked_indi = picked_indi
        self.picked_indj = picked_indj
        self.Num_OperSuccess = Num_OperSuccess
        self.kT = kT
        self.OffSpring = 0

    def ExchangeAtom(self,Parent):
        OffSpring = self.OffSpring
        total_atoms = Parent.get_number_of_atoms()
        times = random.randint(1,total_atoms*3)
        if total_atoms == 1:
            with open('../log.muse','a') as logfile: print >>logfile, "Only 1 atom, no exchange."

            exit(1)
        while OffSpring != Parent:
            OffSpring = copy.deepcopy(Parent)
            while times:
                while True:
                    i,j = random.randint(0,total_atoms-1),random.randint(0,total_atoms-1)
                    if list(OffSpring.get_atomic_numbers())[i] != list(OffSpring.get_atomic_numbers())[j]:
                        k = list(OffSpring.get_atomic_numbers())[i]
                        OffSpring.numbers[i] = list(OffSpring.get_atomic_numbers())[j]
                        OffSpring.numbers[j] = k
                        Read_Write.write_vasp('POSCAR',OffSpring,direct=True,sort=True,vasp5=True)
                        break
                times -=1

        return OffSpring

    def Run_Exch(self):
        spg0 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile: print >>logfile, "Exchange: (%02d: %s)"%(self.picked_indi+1,spg0[0])

            self.ExchangeAtom(self.Old_cry[self.picked_indi][1])
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            #spg1 = Findspg.Findspg(BigDict[nu][ng][enth])

            if enth < self.Old_cry[self.picked_indi][0]:# or spg1[0] != spg0[0]:
                self.Num_OperSuccess[self.ng][7] +=1
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

    def Run_ExchMut(self):
        spg0 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])

        GaussStdDev = float(indict['InitGaussStdDev'][0])
        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Exchange + Mutate %.3f: %02d: %s"%(GaussStdDev,self.picked_indi+1,spg0[0])

            OffSpring = self.ExchangeAtom(self.Old_cry[self.picked_indi][1])
            mutgauss = Mut_Gauss.MUTGAUSS(self.nu,self.ng,self.nn,self.Old_cry,self.BigDict,self.picked_indi,self.V_best,self.Num_OperSuccess,self.kT)
            mutgauss.Mut_Gauss(OffSpring)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)
            with open('../log.muse','a') as logfile: print >>logfile, self.Old_cry[self.picked_indi][0],enth

            #spg1 = Findspg.Findspg(BigDict[nu][ng][enth])
            if enth < self.Old_cry[self.picked_indi][0]:# or spg1[0] != spg0[0]:
                self.Num_OperSuccess[self.ng][8] +=1
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

    def Run_CrossExch(self):
        spg01 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])
        spg02 = Findspg.Findspg(self.Old_cry[self.picked_indj][1])

        #GaussStdDev = float(indict['InitGaussStdDev'][0])
        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Crossover + Exchange: (%02d,%02d: %s, %s)"%(self.picked_indi+1,self.picked_indj+1,spg01[0],spg02[0])

            #crossij = CrossOver.CrossOver(Old_cry[picked_indi][1],Old_cry[picked_indj][1],nu)
            #OffSpring = self.ExchangeAtom(crossij)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile:
                print >>logfile, self.Old_cry[self.picked_indi][0],self.Old_cry[self.picked_indj][0],enth,spg1[0]

            if enth < self.Old_cry[self.picked_indi][0] and enth < self.Old_cry[self.picked_indj][0]:
                self.Num_OperSuccess[self.ng][9] +=1
                ifdo = 0

            else:
                ifdo -= 1
                if indict['Algorithm'][0] == 'MAC' or indict['Algorithm'][0] == 'SA':
                    if np.exp((min(self.Old_cry[self.picked_indi][0],self.Old_cry[self.picked_indj][0]) - enth) / self.kT) < np.random.uniform():
                        if self.Old_cry[self.picked_indi][0] < self.Old_cry[self.picked_indj][0]:
                            BigDict[self.nu][self.ng][self.Old_cry[self.picked_indi][0]] = self.Old_cry[self.picked_indi][1]
                        else:
                            BigDict[self.nu][self.ng][self.Old_cry[self.picked_indj][0]] = self.Old_cry[self.picked_indj][1]

                        #if enth in BigDict[nu][ng].keys():
                            #del BigDict[nu][ng][enth]

            if 111111.11111111 in BigDict[self.nu][self.ng].keys(): del BigDict[self.nu][self.ng][111111.11111111]

        with open('.Num_OperSuccess0','w') as Num_OperSuccessFile:
            pickle.dump(self.Num_OperSuccess,Num_OperSuccessFile)

        os.system("mv .Num_OperSuccess0 .Num_OperSuccess")

        return nk,self.Num_OperSuccess,BigDict
