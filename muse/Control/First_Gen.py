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

import random
from muse.Readwrite import Read_Write
from muse.Calculators import Submit
from muse.Genstrs import GenRandCrystal
from muse.Symmetry import SymmConstr
from muse.Crystaloper import DelRepAndCompl
from muse.Readwrite.ReadInput import indict

def First_Gen(BigDict,nu,ng,IfRep,Maxnu,RepTimes,Old_cry=[]):
    if IfRep:
        Ratio = 0.8
    else:
        Ratio = 1.0

    if nu == int(indict['MinNum_PrimitiveCell'][0]):
        sg_dict = [0]
        for nn0 in range(int(indict['PopSize'][0])):
            nn = nn0 + 1
            if int(indict['IfSymmConstr'][0]) == 1:

                sg_dict,ifempty = SymmConstr.SymmConstr(nu,ng,nn,sg_dict)
                #print sg_dict[1:],len(sg_dict)-1

            elif int(indict['IfSymmConstr'][0]) == 0:
                GenRandCrystal.GenRandCrystal(nu,ng,nn)

            nk,enth,BigDict = Submit.Submit(BigDict,nu,ng,nn,Old_cry)

    else:
        sg_dict = [0]
        for nn0 in range(int(Ratio*int(indict['PopSize'][0]))):
            nn = nn0 + 1
            if int(indict['IfSymmConstr'][0]) == 1:
                sg_dict,ifempty = SymmConstr.SymmConstr(nu,ng,nn,sg_dict)

            elif int(indict['IfSymmConstr'][0]) == 0:
                GenRandCrystal.GenRandCrystal(nu,ng,nn)

            nk,enth,BigDict = Submit.Submit(BigDict,nu,ng,nn,Old_cry)

        if IfRep:
            copieddict = BigDict.copy()

            nf = len(copieddict[Maxnu].keys())
            i=1
            while i < nf:
                copieddict[Maxnu][nf].update(copieddict[Maxnu][i])
                i+=1

            Old_cry = DelRepAndCompl.DelRepeat(copieddict,Maxnu,nf)

            for nn0 in range(int(indict['PopSize'][0])-int(Ratio*int(indict['PopSize'][0]))):

                supercell = Old_cry[nn0][1].repeat(random.choice([(RepTimes,1,1),(1,RepTimes,1),(1,1,RepTimes)]))

                with open('../log.muse','a') as logfile: print >>logfile, "Multiplied..."

                wlabel = indict['NameSys'][0]+": "+str(ng)+'-'+str(nn+1)+" Multiplied."
                Read_Write.write_vasp('POSCAR',supercell,label=wlabel,direct=True,sort=True,vasp5=True)
                nk,enth,BigDict = Submit.Submit(BigDict,nu,ng,nn+1,Old_cry)
                nn +=1

    return nk,BigDict
