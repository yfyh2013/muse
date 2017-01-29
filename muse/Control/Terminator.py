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

from muse.Symmetry import Findspg
from muse.Readwrite.ReadInput import indict

def Terminator(nu,ng,Num_GenBest,Old_cry,BigDict):
    num_best_enth=0
    i=ng-1
    while i >= ng-Num_GenBest+1:
        Old_cry_ng = BigDict[nu][i].items()
        Old_cry_ng.sort()
        
        spg1 = Findspg.Findspg(Old_cry_ng[0][1])
        spg0 = Findspg.Findspg(Old_cry[0][1])
        
        if abs(Old_cry[0][0]-Old_cry_ng[0][0]) < 0.001*Old_cry[0][1].get_number_of_atoms() and spg1[0]==spg0[0]:
            num_best_enth +=1
            
        i -=1
        
    nk = ng
    
    if num_best_enth >= Num_GenBest-1:
        #os.system("cat ../all-enthalpy | grep '-' | sort -k 4 -g >horder")
        with open('../log.muse','a') as logfile: print >>logfile, "Convergence criterion reached."
        nk = int(indict['MaxGen'][0])
        
    return nk

