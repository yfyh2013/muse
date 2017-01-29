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

import os
from muse.Readwrite import Read_Write
from muse.Symmetry import Findspg
from muse.Calculators import Submit
from muse.Readwrite.ReadInput import indict

def DirectOpt(BigDict,Old_cry,nu,ng):
    with open('../log.muse','a') as logfile: print >>logfile
    all_enthfile = open('../all-enthalpy-'+str(nu),'a')
    if int(indict['Num_Keep'][0]) > 0:
        i = 0
        nn = 1
        nkept = 1
        spglist = []
        while nkept <= int(indict['Num_Keep'][0]):
            if int(indict['IfReOptKept'][0]):
                with open('../log.muse','a') as logfile: print >>logfile, "Direct reopt. ..."
                spgnum = Findspg.Findspg(Old_cry[i][1])
                if spgnum[0] not in spglist:
                    spglist.append(spgnum[0])
                    Read_Write.write_vasp('POSCAR',Old_cry[i][1],label=indict['NameSys'][0]+": "+str(ng)+'-'+str(nn),direct=True,sort=True,vasp5=True)

                    nk,enth,BigDict = Submit.Submit(BigDict,nu,ng,nn,Old_cry)

                    nn += 1
                    nkept +=1

            else:
                spgnum = Findspg.Findspg(Old_cry[i][1])
                if spgnum[0] not in spglist:
                    with open('../log.muse','a') as logfile: print >>logfile, "-"*23,"%d-%d"%(ng,nn),"-"*23
                    spglist.append(spgnum[0])
                    with open('../log.muse','a') as logfile:
                        print >>logfile, "%02d: %s, %s %10.4f kept, not reopt."%(i+1,spgnum[0],spgnum[1],Old_cry[i][0])
                        print >>logfile

                    BigDict[nu][ng][Old_cry[i][0]] = Old_cry[i][1].copy()

                    ifexist = os.system("grep %02d-%02d: %s"%(ng,nn,"../all-enthalpy-"+str(nu)))
                    if ifexist != 0:
                        all_enthfile.write(" %02d-%02d:%11s%9s%14.6f%14.6f%14s"%(ng,nn,spgnum[0],spgnum[1],Old_cry[i][0],Old_cry[i][1].get_volume(),'----')+'\n')

                        Read_Write.write_vasp('POSCAR',Old_cry[i][1],label=indict['NameSys'][0]+": "+"%02d-%02d"%(ng,nn)+'  '+spgnum[0]+'  '+str(spgnum[1])+'  '+str(Old_cry[i][0]),direct=True,sort=True,vasp5=True)

                        os.system("cat POSCAR >> ../poscars-%d"%nu)

                    nn += 1
                    nkept +=1

            i +=1

    all_enthfile.close()
    return BigDict
