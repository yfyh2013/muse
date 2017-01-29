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

import sys
from muse.Similarity import Evaluate_Similarity
from muse.Control import Selector
from muse.Symmetry import SymmConstr, Findspg
from muse.Genstrs import GenRandCrystal
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict
from spglib import spglib
from muse.Genstrs import Generate_Supercell


def DelRepeat(BigDict,nu,ng):
    Old_cry,V_best = Selector.GenOrder(BigDict,nu,ng)
    
    spgbase = []
    for jj in range(len(Old_cry)):
        spgnum = Findspg.Findspg(Old_cry[jj][1])
        pc = spglib.find_primitive(Old_cry[jj][1],symprec=1e-1, angle_tolerance=-1.0)
        if pc[0] is not None:
#            print "reduced"
            nat = len(pc[1])
        else:
#            print "not reduced"
            nat = len(Old_cry[jj][1].get_positions())
            
        if [spgnum[0],nat] not in spgbase:
            spgbase.append([spgnum[0],nat])
            
        else:
            spgbase.append(jj)
            ind = spgbase.index([spgnum[0],nat])
            pc1 = spglib.find_primitive(Old_cry[jj][1],symprec=1e-1, angle_tolerance=-1.0)
            pc2 = spglib.find_primitive(Old_cry[ind][1],symprec=1e-1, angle_tolerance=-1.0)
            
            if pc1[0] is not None:
                nat1 = len(pc1[1])
            else:
                nat1 = len(Old_cry[jj][1].get_positions())
            if pc2[0] is not None:
                nat2 = len(pc2[1])
            else:
                nat2 = len(Old_cry[ind][1].get_positions())
            
#            print "nat1= %d, nat2= %d"%(nat1,nat2)
            if pc1[0] is not None and pc2[0] is not None:
                pcell1 = []
                pcell2 = []
                for i in range(len(pc1[0])):
                    pcell1.append([pc1[0][i][j] for j in range(3)])

                for i in range(len(pc2[0])):
                    pcell2.append([pc2[0][i][j] for j in range(3)])
                    
                pos1 = []
                for i in range(len(pc1[1])):
                    pos1.append([pc1[1][i][j] for j in range(3)])

                pos2 = []
                for i in range(len(pc2[1])):
                    pos2.append([pc2[1][i][j] for j in range(3)])
                    
                sys_name1 = [Generate_Supercell.ChemicalSymbols[ii] for ii in pc1[2]]
                sys_name2 = [Generate_Supercell.ChemicalSymbols[ii] for ii in pc2[2]]
                #print sys_name
                primcell1 = Generate_Supercell.Generate_Supercell(sys_name1,cell=pcell1,scaled_positions=pos1)
                primcell2 = Generate_Supercell.Generate_Supercell(sys_name2,cell=pcell2,scaled_positions=pos2)
                
            else:
                primcell1 = Old_cry[jj][1]
                primcell2 = Old_cry[ind][1]
                
            IfTriang = False
            distacc = 0.1
            angacc = 5.0
            r_cutoff = 3.0
            if nat1 == nat2:
                IfTriang = Evaluate_Similarity.evaluate_similarity(primcell1, primcell2, distacc, angacc, r_cutoff)
                
#            spgnum1 = Findspg.Findspg(primcell1)
#            spgnum2 = Findspg.Findspg(primcell2)
#            print spgnum1[0],spgnum2[0]
#            if IfTriang:
#                print "Triangles match."
#            else:
#                print "Triangles mismatch."
            
#            print Old_cry[jj][0] - Old_cry[ind][0]
            if IfTriang:
                del BigDict[nu][ng][Old_cry[jj][0]]
#                print "deleted..."
                
        #print spgbase,spgnum[0],nat
        
    #norepeat = list(set(spgbase))
    #for kk in range(len(norepeat)):
        #count = 0
        #for jj in range(len(Old_cry)):
            #spgnum = Findspg.Findspg(Old_cry[jj][1])
            #if norepeat[kk] == spgnum[0]:
                #count +=1
                #if count >1: del BigDict[nu][ng][Old_cry[jj][0]]#; print "del...%s"%norepeat[kk]
                
    Old_cry=BigDict[nu][ng].items()
    Old_cry.sort()

    sys.stdout.flush()
    #for jj in range(len(Old_cry)):
        #spgnum = Findspg.Findspg(Old_cry[jj][1])
        #print Old_cry[jj][0],spgnum[0],spgnum[1]
        
    return Old_cry


def Complement(BigDict,nu,ng,nn):
    Old_cry = DelRepeat(BigDict,nu,ng)
    if len(Old_cry) < int(indict['PopSize'][0]):
        with open('../log.muse','a') as logfile: print >>logfile, "Complement some random strs..."
        
        nn = len(Old_cry)+1
        sg_dict = [0]
        while nn <= int(indict['PopSize'][0]):
            if int(indict['IfSymmConstr'][0]) == 1:
                SymmConstr.SymmConstr(nu,ng,nn,sg_dict)

            elif int(indict['IfSymmConstr'][0]) == 0:
                GenRandCrystal.GenRandCrystal(nu,ng,nn)

            nk,enth,BigDict = Submit.Submit(BigDict,nu,ng,nn+int(indict['PopSize'][0])-len(Old_cry),Old_cry)

            nn +=1
            
        Old_cry=BigDict[nu][ng].items()
        Old_cry.sort()

    return Old_cry
