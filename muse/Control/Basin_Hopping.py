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
import random,copy,Selector
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Calculators import DirectOpt,Submit
from muse.Readwrite.ReadInput import indict

def Disp_Atoms(atm,dr=0.1):
    atoms = copy.deepcopy(atm)
    pos_old = atoms.get_positions()
    disp = np.random.uniform(-1, 1., (len(atoms), 3))
    pos_new = pos_old + dr * disp
    atoms.set_positions(pos_new)
    return atoms

def Basin_Hopping(BigDict,nu,ng,kT):
    nb = 0
    kT = float(indict['InitTemp_anneal'][0]) * 1.380658e-23 / 1.60217733e-19
    
    while nb < int(indict['Num_GenBest'][0]):
        kT *= 0.90
        
        BigDict[nu][ng+nb] = {}
        O_cry,V_best = Selector.GenOrder(BigDict,nu,ng+nb-1)
        Old_cry = copy.deepcopy(O_cry)
        BigDict = DirectOpt.DirectOpt(BigDict,copy.deepcopy(Old_cry),nu,ng+nb)
        
        for nn in range(len(O_cry)):
            if nn == int(indict['PopSize'][0])-int(indict['Num_Keep'][0])-1: break
            #picked_indi,picked_indj = Selector.IndiSelector(BigDict,nu,ng+nb-1)
            spg0 = Findspg.Findspg(Old_cry[nn][1])
            with open('../log.muse','a') as logfile: print >>logfile, "Basin hopping for %02d: %s"%(nn+1,spg0[0])
            
            ifdo = 1
            while ifdo:
                dispatoms = Disp_Atoms(Old_cry[nn][1],dr=0.5*(1+random.random()))
                wlabel = indict['NameSys'][0]+": "+str(ng+nb)+'-'+str(nn+int(indict['Num_Keep'][0])+1)
                Read_Write.write_vasp('POSCAR',dispatoms,label=wlabel,direct=True,sort=True,vasp5=True)
                nk,enth,BigDict = Submit.Submit(BigDict,nu,ng+nb,nn+int(indict['Num_Keep'][0])+1,Old_cry)
                
                spg1 = Findspg.Findspg(BigDict[nu][ng+nb][enth])
                
                with open('../log.muse','a') as logfile: print >>logfile, Old_cry[nn][0],enth,spg1[0]
                
                if enth < Old_cry[nn][0]:# or spg1[0] != spg0[0]:
                    ifdo = 0
                    
                else:
                    ifdo -= 1
                    
                    if np.exp((Old_cry[nn][0] - enth) / kT) < np.random.uniform():
                        BigDict[nu][ng][Old_cry[nn][0]] = Old_cry[nn][1]
                        
                        #if enth in BigDict[nu][ng].keys():
                        #del BigDict[nu][ng][enth]
                        
                if 111111.11111111 in BigDict[nu][ng].keys(): del BigDict[nu][ng][111111.11111111]
                
        nb +=1
        
    return BigDict

if __name__ == '__main__':
    Basin_Hopping()
