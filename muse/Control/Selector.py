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

import random,copy
from muse.Readwrite.ReadInput import indict

if len(indict['AtomicNumber']) > 1:
    opitems = ['Mut','Cross','CrossMut','Slip','RandomMove','Ripple','Twist','Exch','ExchMut','CrossExch']
        
else:
    opitems = ['Mut','Cross','CrossMut','Slip','RandomMove','Ripple','Twist']

    
initprob = [1.0/len(opitems) for _ in opitems]

def Random_Select(items,prob):
    x = random.uniform(0,1)
    cp = 0.0
    
    for item, ip in zip(items,prob):
        cp += ip
        if x < cp: break
    
    return item
    
def GenOrder(BigDict,nu,ng):
    Dict = copy.deepcopy(BigDict)
    O_cry=Dict[nu][ng].items()
    O_cry.sort()
    
    Old_cry = copy.deepcopy(O_cry)
    V_best = Old_cry[0][1].get_volume()
    
    return Old_cry,V_best
    
def SelectProb(BigDict,nu,ng):
    (Old_cry,V_best) = GenOrder(BigDict,nu,ng)
    indprob = []
    finalindprob = []
    for m in range(len(Old_cry)):
        if Old_cry[len(Old_cry)-1][0]-Old_cry[0][0] == 0:
            with open('../log.muse','a') as logfile:
                print >>logfile, "The enthalpy differences are so small, maybe converged otherwise increase PopSize."
            
            exit(0)
            
        else:
            select_prob = 1-(Old_cry[m][0]-Old_cry[0][0])/(Old_cry[len(Old_cry)-1][0]-Old_cry[0][0])
            indprob.append(select_prob)
            
    finalindprob = [round(i / sum(indprob),3) for i in indprob]
    
    return finalindprob
    
def OperSelector(prob,nu,ng):
    inprob = initprob
    
    if ng >=3:
        inprob = prob
        
    picked_operator = Random_Select(opitems,inprob)
    
    return picked_operator
    
def IndiSelector(BigDict,nu,ng):
    
    inditems = [i for i in range(int(round(float(indict['Perc4NextGen'][0])*int(indict['PopSize'][0]))))]
    
    finalindprob = SelectProb(BigDict,nu,ng)
    picked_indi = Random_Select(inditems,finalindprob)
    picked_indj = Random_Select(inditems,finalindprob)
    
    max_try = 0
    while picked_indi == picked_indj:
        max_try += 1
        picked_indj = Random_Select(inditems,finalindprob)
        if max_try >100:
            break
        
    Old_cry,V_best = GenOrder(BigDict,nu,ng)
    
    return picked_indi,picked_indj
