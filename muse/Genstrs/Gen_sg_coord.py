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
from muse.Symmetry import Wyckoff_dict
from muse.Genstrs import FindRatio

def IfEnoughFar(pos1,pos2):
    allm = []
    for i in range(3):
        allm.append(abs(pos1[i]-pos2[i]))
        
    return max(allm) > 0.1

def Gen_sg_coord(sg_num,Num_atomsInChemForm):
    mm = []
    colist = []
    for i in range(len(Wyckoff_dict.Wyck(np.random.uniform(0.1, 0.9, (19, 3)))[sg_num].items())):
        mm.append(Wyckoff_dict.Wyck(np.random.uniform(0.1, 0.9, (19, 3)))[sg_num].items()[i][1][0])
    
    mm.sort()
    #print sg_num
    #print mm
    m = 1
    issuccess = False
    colist =[]
    Nspec = []
    while sum(Nspec) <= sum(mm) and not issuccess:
        #print "m = %d"%m
        Nspec = [ i*m for i in Num_atomsInChemForm]
#        print "Nspec:",Nspec
        
        occ,aa,ratiolist,issuccess = FindRatio.FindRatio(Nspec,mm)
#        print "issuccess:",issuccess
        #print Nspec,Num_atomsInChemForm
        #print ratiolist
        if issuccess:
            ifok = [0]
            
            while 0 in ifok:
                #print "="*20
                ifok = []
                wdict = Wyckoff_dict.Wyck(np.random.uniform(0.1, 0.9, (19, 3)))[sg_num]
                allpos = []
                for k in wdict.keys():
                    allpos.append(wdict.values()[k-1][1])
                    
                for pi,pos1 in enumerate(allpos):
                    for pj,pos2 in enumerate(allpos):
                        if pi < pj:
                            if not IfEnoughFar(pos1,pos2):
                                ifok.append(0)
                #print "ifok:", ifok    
            for i in ratiolist:
                colist.append(wdict[len(mm)-i][1])
    
        m +=1
#    print "colist:",colist
    return occ,colist,aa,issuccess

if __name__ == '__main__':

    Num_atomsInChemForm = [1,1,3]
    occ,colist,Nspec,issuccess = Gen_sg_coord(63,Num_atomsInChemForm)
    print "testing..."
    print occ,colist,Nspec,issuccess
