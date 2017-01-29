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

import os,random,copy
import numpy as np
from spglib import spglib
from muse.Readwrite.ReadInput import indict
from muse.Symmetry.Const import *
from muse.Symmetry.SG import *
from muse.Genstrs import GENCELL,Gen_sg_coord
from muse.Genstrs.GenPC import GenPC
from muse.Genstrs import Generate_Supercell

def Run_genpc(nu,sg_dict):
    Num_atomsInChemForm = [int(i) for i in indict['Num_atomsInChemForm']]
    oldlist = [nu*i for i in Num_atomsInChemForm]

    n = 100
    iffind = False
    Isdirectmatch = False
    while not iffind:
        issuccess = False
        colist = []
        while not issuccess or colist == []:
            if os.path.isfile('../.sg'):
                with open('../.sg','r') as sgfile:
                    n = int(sgfile.readline()[0])

            if 'KnownSGN' in indict and n < len(indict['KnownSGN']):
                sg_num = int(indict['KnownSGN'][n])
                with open('../.sg','w') as sgfile:
                    print >>sgfile, n+1
            else:
                sg_num = random.randint(2,230)

            #print "sg:",sg_num
            occ,colist,Nspec,issuccess = Gen_sg_coord.Gen_sg_coord(sg_num,Num_atomsInChemForm)
            
        
#        print sg_num,occ,colist
        gencell2 = GENCELL.GENCELL(sg_num)
        cellp = gencell2.Gen_sg_latt()
        
        while min(cellp[0],cellp[1],cellp[2]) < 1.0:
            for i in range(3): cellp[i] *= 10
        
        #occ = [{indict['NameOfAtoms'][i]:Nspec[i]} for i in range(len(indict['NameOfAtoms']))]
        ineqsites = [list(colist[i]) for i in range(len(colist))]
        genpc = GenPC(sg_num,cellp,occ,ineqsites)
        
        ## Not reduce:
        primcell = genpc.getCrystalStructure(IfReduce=False)
        newlist = []
        for a in [i for i in indict['NameOfAtoms']]:
            newlist.append(primcell.get_chemical_symbols().count(a))
        
        ifone = []
        ifdive = []
        for i,num in enumerate(oldlist):
            ifone.append(float(oldlist[i]/newlist[i]))
            ifdive.append(oldlist[i]%newlist[i] == 0)
            
        if len(set(ifone)) == 1.0 and ifone[0] >= 1.0 and len(set(ifdive)) == 1.0 and ifdive[0]:
#            print "Succeeded after not reduced."
            iffind = True
            if sum(oldlist) == sum(newlist):
                Isdirectmatch = True
        else:
#            print "Failed after not reduced."
#            print oldlist
#            print newlist
            iffind = False
        
        ## Reduce:
        if not iffind:
            #print "Try reduce:"
            primcellbak = copy.deepcopy(primcell)
            primcell = spglib.find_primitive(primcell,symprec=0.01,angle_tolerance=0.01)
            if primcell is not None:
                pcell = []
                for i in range(len(primcell[0])):
                    pcell.append([primcell[0][i][j] for j in range(3)])
                
                V_init = abs(np.linalg.det(pcell))
                for ll in range(len(pcell)):
                    for jj in range(3):
                        pcell[ll][jj] *= (nu*float(indict['V_guess'][0])/V_init)**(1./3.)
                        
                pos = []
                for i in range(len(primcell[1])):
                    pos.append([primcell[1][i][j] for j in range(3)])
                
                sys_name = [Generate_Supercell.ChemicalSymbols[ii] for ii in primcell[2]]
                #print sys_name
                primcell = Generate_Supercell.Generate_Supercell(sys_name,cell=pcell,scaled_positions=pos)
                
                newlist = []
                for a in [i for i in indict['NameOfAtoms']]:
                    newlist.append(primcell.get_chemical_symbols().count(a))
                
                ifone = []
                ifdive = []
                for i,num in enumerate(oldlist):
                    ifone.append(float(oldlist[i]/newlist[i]))
                    ifdive.append(oldlist[i]%newlist[i] == 0)
                    
                if len(set(ifone)) == 1.0 and ifone[0] >= 1.0 and len(set(ifdive)) == 1.0 and ifdive[0]:
#                    print "Succeeded after reduced."
#                    print oldlist
#                    print newlist
                    iffind = True
                else:
#                    print "Failed after reduced."
#                    print oldlist
#                    print newlist
                    iffind = False
            else:
                primcell = primcellbak
    return primcell,sg_num,Isdirectmatch

if __name__ == '__main__':
    pc = Run_genpc(3,[0])
    print pc
