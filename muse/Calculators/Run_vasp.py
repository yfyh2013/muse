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
import os,subprocess,time,datetime
from muse.Readwrite import Read_Write,MakeInputfile
from muse.Readwrite.ReadInput import indict


def Run_vasp(nu,ng,nn,ParSubCom):
    Kfile = open('../kgrids-'+str(nu),'a')
    Nopt = int(os.popen("grep Opt-step ../../INCARs | wc -l").readline())
    
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                 int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
                 
    atomorder.sort()
    os.system("cp ../../POTCAR-%s POTCAR"%atomorder[0][0])
    
    i = 1
    while i < len(atomorder):
        os.system("cat ../../POTCAR-%s >> POTCAR"%atomorder[i][0])
        i +=1
        
    #os.system("cp ../POTCAR .")
    i = 1
    while i <= Nopt:
        K=[]
        New_cry = Read_Write.read_vasp('POSCAR')
        cl = New_cry.get_cell()
        
        # Produce the K Grid for the str.
        for j in range(3):
            K.append(1+int(1.0/(float(indict['K_res'][i-1])*(cl[j][0]**2+cl[j][1]**2+cl[j][2]**2)**.5)))
            
        Kfile.write("%02d-%02d-%02d: %2d %2d %2d"%(ng,nn,i,K[0],K[1],K[2])+'\n')
        Kfile.flush()
        
        # Write files and running.
        MakeInputfile.MakeInputfile('../../INCARs','INCAR','Opt-step',i)
        
        os.system("sed -i 's/From-muse.in/%.2f/g' INCAR"%(float(indict['Pressure'][0])*10))
        os.system("echo A > KPOINTS")
        os.system("echo 0 >> KPOINTS")
        os.system("echo Gamma >> KPOINTS")
        os.system("echo %d  %d  %d >> KPOINTS"%(K[0],K[1],K[2]))
        os.system("echo 0  0  0 >> KPOINTS")
        #os.system("%s"%(ParSubCom))
        
        start = datetime.datetime.now()
        go = subprocess.Popen(ParSubCom,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        while go.poll() is None:
            time.sleep(5)
            now = datetime.datetime.now()
            if (now - start).seconds > float(indict['MaxHour'][0])*3600:
                os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
                break

        if int(os.popen("cat CONTCAR | wc -l").readline()) > 0:
            os.system("cp CONTCAR POSCAR")
        with open('../log.muse','a') as logfile: print >>logfile, "Opt. %d ... done."%i
        
        i += 1
        
    try:
        enth = float(os.popen("grep enthalpy OUTCAR | tail -1 | awk '{print $5}'").readline())
        enth = float("%.8f"%enth)
        FinalPressure = float(os.popen("grep external OUTCAR | awk '{print $9/10}'").readline())
        os.system("rm OUTCAR")
        
    except:
        enth = 111111.11111111
        FinalPressure = 111111.11111111
        pass
    
    return enth,FinalPressure
