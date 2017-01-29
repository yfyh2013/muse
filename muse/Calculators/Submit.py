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

import os,copy,pickle
from muse.Symmetry import Findspg
from muse.Crystaloper import HardConstrain
from muse.Control import Terminator
from muse.Calculators import Run_vasp,Run_siesta,Run_pwscf,Run_abinit,Run_castep,Run_lammps
from muse.Readwrite import Read_Write
from muse.Readwrite.ReadInput import indict

#Parallel submiting command.
ParSubCom = ''
for i in range(len(indict['ParSubCom'])):
    ParSubCom += indict['ParSubCom'][i]
    ParSubCom += ' '

def Submit(BigDict,nu,ng,nn,Old_cry):

    # Open files to write.
    all_enthfile = open('../all-enthalpy-'+str(nu),'a')
    best_enthfile = open('../best-enthalpy-'+str(nu),'a')
    RestartFile = open('../.Restart0','w')

    with open('../log.muse','a') as logfile:
        if nu < 10:
            print >>logfile, str(nu)*22,"("+"%d-%d"%(ng,nn)+")",str(nu)*22
        else:
            print >>logfile, str(nu)*11,"("+"%d-%d"%(ng,nn)+")",str(nu)*11

        print >>logfile, "Local opt. ..."

    if int(indict['OptCode'][0]) == 1:
        enth,FinalPressure = Run_vasp.Run_vasp(nu,ng,nn,ParSubCom)

    elif int(indict['OptCode'][0]) == 2:
        enth,FinalPressure = Run_siesta.Run_siesta(nu,ng,nn,ParSubCom)

    elif int(indict['OptCode'][0]) == 3:
        enth,FinalPressure = Run_pwscf.Run_pwscf(nu,ng,nn,ParSubCom)

    elif int(indict['OptCode'][0]) == 4:
        enth,FinalPressure = Run_castep.Run_castep(nu,ng,nn,ParSubCom)

#    elif int(indict['OptCode'][0]) == 5:
#        enth,FinalPressure = Run_abinit.Run_abinit(nu,ng,nn,ParSubCom)

    elif int(indict['OptCode'][0]) == 5:
        enth,FinalPressure = Run_lammps.Run_lammps(nu,ng,nn,ParSubCom)

    if ng >= 2 and nn == 1 and int(indict['IfCompeteBetweenGens'][0]):
        BigDict[nu][ng] = copy.deepcopy(BigDict[nu][ng-1])

    if int(os.popen("cat CONTCAR | wc -l").readline()) > 0:
        opted_POS = Read_Write.read_vasp('CONTCAR')
        L_dangle, L_angle, L_len, L_dist = HardConstrain.HardConstrain(opted_POS, nu, 35, 145)

    else:
        L_dangle = True

    if not L_dangle and not L_angle and not L_len and not L_dist and enth != 111111.11111111:
        spgnum = Findspg.Findspg(opted_POS)
        Read_Write.write_vasp('POSCAR',opted_POS,label=indict['NameSys'][0]+": "+"%02d-%02d"%(ng,nn)+'  '\
                               +spgnum[0]+'  '+str(spgnum[1])+'  '+str(enth), direct=True,sort=True,vasp5=True)

        os.system("cat POSCAR >> ../poscars-%d"%nu)
        BigDict[nu][ng][enth]=Read_Write.read_vasp('CONTCAR')
        spgnum = Findspg.Findspg(BigDict[nu][ng][enth])

        all_enthfile.write(" %02d-%02d:%11s%9s%14.6f%14.6f%14.6f"%
                          (ng,nn,spgnum[0],spgnum[1],
                           enth,BigDict[nu][ng][enth].get_volume(),FinalPressure)+'\n')

        all_enthfile.flush()

        pickle.dump(BigDict,RestartFile)
        RestartFile.close()
        os.system("mv ../.Restart0 ../.Restart")
        with open('../log.muse','a') as logfile:
            print >>logfile, "%11s%9s%14.6f%14.6f"%(spgnum[0],spgnum[1],enth,BigDict[nu][ng][enth].get_volume())
            print >>logfile

    else:
        with open('../log.muse','a') as logfile: print >>logfile, "Nonphysical str. discarded."
        enth = 111111.11111111
        BigDict[nu][ng][enth] = Read_Write.read_vasp('POSCAR')

    nk = ng
    if nn == int(indict['PopSize'][0]):
        #all_enthfile.write("="*58+'\n')
        Old_cry=BigDict[nu][ng].items()
        Old_cry.sort()
        spgnum = Findspg.Findspg(BigDict[nu][ng][Old_cry[0][0]])
        ifexist = os.system("grep %02d: %s"%(ng,"../best-enthalpy-"+str(nu)))
        if ifexist != 0:
            best_enthfile.write("%02d: %11s%9s%14.6f%14.6f%14.6f"%
                               (ng,spgnum[0],spgnum[1],Old_cry[0][0],Old_cry[0][1].get_volume(),FinalPressure)+'\n')

            best_enthfile.flush()

        V_best = Old_cry[0][1].get_volume()
        with open('../log.muse','a') as logfile: print >>logfile, "V_best: %.4f"%V_best
        #print Old_cry
        if ng >= int(indict['Num_GenBest'][0]):
            nk = Terminator.Terminator(nu,ng,int(indict['Num_GenBest'][0]),Old_cry,BigDict)

    return nk,enth,BigDict
