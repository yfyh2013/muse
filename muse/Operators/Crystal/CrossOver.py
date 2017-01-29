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
from muse.Crystaloper import Centering,HardConstrain
from muse.Symmetry import Findspg
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict

class CROSSOVER:
    def __init__(self,nu,ng,nn,Old_cry,BigDict,picked_indi,picked_indj,V_best,Num_OperSuccess,kT):
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

    def Exchange_Parents(self,Parent1,Parent2):
        ParentTmp  = copy.deepcopy(Parent1)
        Parent1  = copy.deepcopy(Parent2)
        Parent2  = copy.deepcopy(ParentTmp)
        return Parent1,Parent2

    def RandomCutandPaste(self,Parent1,Parent2,ChooseParent,OffSpring,randaxis,randpos,randcomp):
        if ChooseParent == 2:
            (Parent1,Parent2) = self.Exchange_Parents(Parent1,Parent2)

        #pos1= Parent1.get_scaled_positions()
        pos2= Parent2.get_scaled_positions()
        pos = OffSpring.get_scaled_positions()

        AtomNum = Parent2.numbers
        choosed = []
        pos2index = []
        for i in range(len(pos2)):
            if randcomp:
                if pos2[i][randaxis] >= randpos:
                    pos2index.append(i)
                    choosed.append(pos2[i])

            if not randcomp:
                if pos2[i][randaxis] <= randpos:
                    pos2index.append(i)
                    choosed.append(pos2[i])

        j=0
        for i in range(len(pos)):
            if randcomp:
                if pos[i][randaxis] >= randpos:
                    if j > len(choosed)-1: break
                    pos[i]=choosed[j]
                    OffSpring.numbers[i]=AtomNum[int(pos2index[j])]
                    j+=1
            if not randcomp:
                if pos[i][randaxis] <= randpos:
                    if j > len(choosed)-1: break
                    pos[i]=choosed[j]
                    OffSpring.numbers[i]=AtomNum[int(pos2index[j])]
                    j+=1

        if ChooseParent == 2:
            (Parent1,Parent2) = self.Exchange_Parents(Parent1,Parent2)

        return pos

    def CrossOver(self,Par1,Par2,nu):
        Parent1 = Centering.Centering(Par1)
        Parent2 = Centering.Centering(Par2)

        total_atoms = Parent1.get_number_of_atoms()
        Max_try = 0
        L_len=True; L_dist=True; L_angle=True; L_dangle=True

        ChooseParent = random.randint(1, 2)
        if ChooseParent == 1:
            OffSpring = copy.deepcopy(Parent1)
        if ChooseParent == 2:
            OffSpring = copy.deepcopy(Parent2)

        while L_dangle or L_len or L_dist or L_angle:
            if total_atoms == 1: break

            newcell=[]
            for i in range(3):
                rand=random.random()
                newcell.append(rand*Parent1.get_cell()[i]+(1.0-rand)*Parent2.get_cell()[i])

            OffSpring.set_cell((newcell),scale_atoms=True)

            randaxis = random.randint(0,2)
            #randpos = random.gauss(0.5,0.1)
            randpos = 0.25 + 0.5*random.random()
            randcomp = random.randint(0,1)

            pos = self.RandomCutandPaste(Parent1,Parent2,ChooseParent,OffSpring,randaxis,randpos,randcomp)

            OffSpring.set_scaled_positions(pos)

            L_dangle,L_angle,L_len,L_dist = HardConstrain.HardConstrain(OffSpring,nu)
            Max_try += 1

            if Max_try >= 100: break

        #pos = OffSpring.get_scaled_positions()
        Read_Write.write_vasp('POSCAR',OffSpring,direct=True,sort=True,vasp5=True)
        sorted_atoms = Read_Write.read_vasp('POSCAR')

        aorder = {indict['NameOfAtoms'][i]:[nu*int(indict['Num_atomsInChemForm'][i]),int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}
        atomorder = aorder.items()

        atomorder.sort()

        y = []
        for i in range(len(atomorder)):
            for _ in range(int(atomorder[i][1][0])):
                y.append(int(atomorder[i][1][1]))

        #atomic_num = []
        #for i in sorted_atoms.get_atomic_numbers():
            #atomic_num.append(i)

        for i in range(len(y)):
            sorted_atoms.numbers[i] = y[i]

        Read_Write.write_vasp('POSCAR',sorted_atoms,direct=True,sort=True,vasp5=True)

        return sorted_atoms

    def Run_Cross(self):
        spg01 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])
        spg02 = Findspg.Findspg(self.Old_cry[self.picked_indj][1])

        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Crossover: (%02d,%02d: %s, %s)"%(self.picked_indi+1,self.picked_indj+1,spg01[0],spg02[0])

            self.CrossOver(self.Old_cry[self.picked_indi][1],self.Old_cry[self.picked_indj][1],self.nu)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile:
                print >>logfile, self.Old_cry[self.picked_indi][0],self.Old_cry[self.picked_indj][0],enth,spg1[0]

            if enth < self.Old_cry[self.picked_indi][0] and enth < self.Old_cry[self.picked_indj][0]:
                self.Num_OperSuccess[self.ng][1] +=1
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

    def Run_CrossMut(self):

        spg01 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])
        spg02 = Findspg.Findspg(self.Old_cry[self.picked_indj][1])

        #GaussStdDev = float(indict['InitGaussStdDev'][0])
        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Crossover + Mutate: (%02d,%02d: %s, %s)"%(self.picked_indi+1,self.picked_indj+1,spg01[0],spg02[0])

            crossij = self.CrossOver(self.Old_cry[self.picked_indi][1],self.Old_cry[self.picked_indj][1],self.nu)
            mutgauss = Mut_Gauss.MUTGAUSS(self.nu,self.ng,self.nn,self.Old_cry,self.BigDict,self.picked_indi,self.V_best,self.Num_OperSuccess,self.kT)
            mutgauss.Mut_Gauss(crossij)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile:
                print >>logfile, self.Old_cry[self.picked_indi][0],self.Old_cry[self.picked_indj][0],enth,spg1[0]

            if enth < self.Old_cry[self.picked_indi][0] and enth < self.Old_cry[self.picked_indj][0]:
                self.Num_OperSuccess[self.ng][2] +=1
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

    def Run_MutCross(self):
        spg01 = Findspg.Findspg(self.Old_cry[self.picked_indi][1])
        spg02 = Findspg.Findspg(self.Old_cry[self.picked_indj][1])

        #GaussStdDev = float(indict['InitGaussStdDev'][0])
        ifdo = int(indict['MaxTry_EachOperator'][0])
        while ifdo:
            with open('../log.muse','a') as logfile:
                print >>logfile, "Mutate + Crossover: (%02d,%02d: %s, %s)"%(self.picked_indi+1,self.picked_indj+1,spg01[0],spg02[0])

            muti = Mut_Gauss.MUTGAUSS(self.nu,self.ng,self.nn,self.Old_cry[self.picked_indi][1],self.BigDict,self.picked_indi,self.V_best,self.Num_OperSuccess,self.kT)
            mutj = Mut_Gauss.MUTGAUSS(self.nu,self.ng,self.nn,self.Old_cry[self.picked_indj][1],self.BigDict,self.picked_indi,self.V_best,self.Num_OperSuccess,self.kT)
            Muti = muti.Mut_Gauss()
            Mutj = mutj.Mut_Gauss()

            self.CrossOver(Muti,Mutj,self.nu)
            nk,enth,BigDict = Submit.Submit(self.BigDict,self.nu,self.ng,self.nn,self.Old_cry)

            spg1 = Findspg.Findspg(BigDict[self.nu][self.ng][enth])

            with open('../log.muse','a') as logfile:
                print >>logfile, self.Old_cry[self.picked_indi][0],self.Old_cry[self.picked_indj][0],enth,spg1[0]

            if enth < self.Old_cry[self.picked_indi][0] and enth < self.Old_cry[self.picked_indj][0]:
                self.Num_OperSuccess[self.ng][3] +=1
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
