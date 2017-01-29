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


def Run_abinit(nu,ng,nn,ParSubCom):
    Nopt = int(os.popen("grep Opt-step ../../abinit.in | wc -l").readline())
    
    atomorder = {indict['NameOfAtoms'][i]:\
                 [nu*int(indict['Num_atomsInChemForm'][i]),\
                 int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
                     
    atomorder.sort()
    
    os.system("cp ../../files .")
    indict['NameOfAtoms'].sort()
    for xx in indict['NameOfAtoms']:
        os.system("echo %s.pot >>files"%xx)
        os.system("cp ../../%s.pot ."%xx)
        
    atoms = Read_Write.read_vasp("POSCAR")
    Read_Write.write_abinit("latt.in",atoms)
    os.system("cp latt.in latt0.in")
    
    i = 1
    while i <= Nopt:
        K=[]
        Kfile = open('../kgrids-'+str(nu),'a')
        New_cry = Read_Write.read_vasp('POSCAR')
        cl = New_cry.get_cell()
        
        # Produce the K Grid for the str.
        for j in range(3):
            K.append(1+int(1.0/(float(indict['K_res'][i-1])*(cl[j][0]**2+cl[j][1]**2+cl[j][2]**2)**.5)))
            
        Kfile.write("%02d-%02d-%02d: %2d %2d %2d"%(ng,nn,i,K[0],K[1],K[2])+'\n')
        Kfile.flush()
        
        with open('latt.in','a') as includefile:
            print >>includefile, "ngkpt      %2d %2d %2d"%(K[0],K[1],K[2])
        
        # Write files and running.
        MakeInputfile.MakeInputfile('../../abinit.in','abinit.in','Opt-step',i)
        
        os.system("cat latt.in >> abinit.in")
        
        os.system("sed -i 's/From-muse.in/%.12f/g' abinit.in"%(-3.39892892272e-05*float(indict['Pressure'][0])))
        
        start = datetime.datetime.now()
        go = subprocess.Popen(ParSubCom,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        while go.poll() is None:
            time.sleep(2)
            now = datetime.datetime.now()
            if (now - start).seconds > float(indict['MaxHour'][0])*3600:
                os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
                break
            
        try:
            New_cry = Read_Write.read_abinit("OUT",nu)
            Read_Write.write_abinit("latt.in",New_cry)
            os.system("cp latt.in latt0.in")
            
        except:
            os.system("cp latt0.in latt.in")
            pass
        
        os.system("cp OUT OUT%s"%i)
        os.system("rm abinit.out*")
        with open('../log.muse','a') as logfile: print >>logfile, "Opt. %d ... done."%int(i)
        
        i += 1
        
    New_cry = Read_Write.read_abinit("OUT",nu)
    Read_Write.write_vasp("CONTCAR",New_cry,label=indict['NameSys'][0]+": "+str(ng)+'-'+str(nn),\
                             direct=True,sort=True,vasp5=True)
                             
    #os.system("cp CONTCAR CONTCAR%s-%s"%(str(ng),str(nn)))
    
    try:
        Etol = float(os.popen("grep etotal OUT | tail -1 | awk '{print $5}'").readline())
        FinalPressure = float(os.popen("grep Pressure OUT | awk '{print $8}'").readline())
        V = float(os.popen("grep \(ucvol\) OUT | tail -1 | awk '{print $6}'").readline())
        enth = Etol + FinalPressure * V / 29421.033
        enth = float("%.8f"%enth*27.21138)
        
    except:
        enth = 111111.11111111
        FinalPressure = 111111.11111111
        pass
    
    return enth,FinalPressure
