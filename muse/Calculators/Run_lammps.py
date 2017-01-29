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
from muse.Readwrite import Read_Write
from muse.Readwrite.ReadInput import indict

def Run_lammps(nu,ng,nn,ParSubCom):
    os.system("cp ../../lammps.in .")
    os.system("cp ../../*.pot . 2>>/dev/null")
    New_cry = Read_Write.read_vasp('POSCAR')
    Read_Write.write_lammps("data",New_cry)
    
    os.system("sed -i 's/From-muse.in/%.2f/g' lammps.in"%(10000*float(indict['Pressure'][0])))
    
    start = datetime.datetime.now()
    go = subprocess.Popen(ParSubCom,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    while go.poll() is None:
        time.sleep(2)
        now = datetime.datetime.now()
        if (now - start).seconds > float(indict['MaxHour'][0])*3600:
            os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
            break
    
    with open('../log.muse','a') as logfile: print >>logfile, "Opt. %d ... done."%int(1)
    
    #os.system("cp OUT OUT%s-%s"%(str(ng),str(nn)))
    try:
        os.system("cp POSCAR CONTCAR")
        New_cry = Read_Write.read_lammps("relax.lammpstrj",nu)
        Read_Write.write_vasp("CONTCAR",New_cry,label=indict['NameSys'][0]+": "+str(ng)+'-'+str(nn),\
                             direct=True,sort=True,vasp5=True)
        
        enth = float(os.popen("grep Enthalpy= OUT | awk '{print $2}'").readline())
        enth = float("%.8f"%enth)
        FinalPressure = float(os.popen("grep Pressure= OUT | awk '{print $2/10000}'").readline())
        #os.system("rm OUT")
        
    except:
        enth = 111111.11111111
        FinalPressure = 111111.11111111
        pass
    
    return enth,FinalPressure
