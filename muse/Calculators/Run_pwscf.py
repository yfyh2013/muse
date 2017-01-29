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


def Run_pwscf(nu,ng,nn,ParSubCom):
    Nopt = int(os.popen("grep Opt-step ../../pwscf.in | wc -l").readline())
    
    atomorder = {indict['NameOfAtoms'][i]:\
                [nu*int(indict['Num_atomsInChemForm'][i]),\
                 int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
                 
    atomorder.sort()
    
    for xx in indict['NameOfAtoms']:
        os.system("cp ../../%s.UPF ."%(xx))
        
    atoms = Read_Write.read_vasp("POSCAR")
    Read_Write.write_pwscf("latt.in",atoms)
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
            print >>includefile, "K_POINTS automatic"
            print >>includefile, "%2d %2d %2d %2d %2d %2d"%(K[0],K[1],K[2],0,0,0)
            
        # Write files and running.
        MakeInputfile.MakeInputfile('../../pwscf.in','pwscf.in','Opt-step',i)
        
        os.system("cat latt.in >> pwscf.in")
        os.system("sed -i 's/PRESS/%.2f/g' pwscf.in"%(10*float(indict['Pressure'][0])))
        os.system("sed -i 's/NAT/%.0f/g' pwscf.in"%(nu*sum([int(x) for x in indict['Num_atomsInChemForm']])))
        os.system("sed -i 's/NTYP/%.0f/g' pwscf.in"%(len(indict['Num_atomsInChemForm'])))
        
        start = datetime.datetime.now()
        go = subprocess.Popen(ParSubCom,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        
        while go.poll() is None:
            time.sleep(2)
            now = datetime.datetime.now()
            if (now - start).seconds > float(indict['MaxHour'][0])*3600:
                os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
                break
            
        try:
            New_cry = Read_Write.read_pwscf("OUT",nu)
            Read_Write.write_pwscf("latt.in",New_cry)
            os.system("cp latt.in latt0.in")
            
        except:
            os.system("cp latt0.in latt.in")
            pass
        
        #os.system("cp OUT OUT%s"%(str(i)))
        with open('../log.muse','a') as logfile: print >>logfile, "Opt. %d ... done."%int(i)
        
        i += 1
        
        New_cry = Read_Write.read_pwscf("OUT",nu)
        Read_Write.write_vasp("CONTCAR",New_cry,label=indict['NameSys'][0]+": "+str(ng)+'-'+str(nn),\
                              direct=True,sort=True,vasp5=True)
        
        #os.system("cp CONTCAR CONTCAR%s-%s"%(str(ng),str(nn)))
        
        try:
            Etol = float(os.popen("grep ! OUT | tail -1 | awk '{print $5}'").readline())
            FinalPressure = float(os.popen("grep P= OUT | tail -1 | awk '{print $6/10}'").readline())
            V = float(os.popen("grep final OUT | grep volume | awk '{print $5}'").readline())
            enth = Etol + FinalPressure * V * 0.529177**3*0.01 / 1.60217653
            enth = float("%.8f"%enth)
            #with open('test.log','w') as testfile: print >>testfile, "Etol= %10.4f, P= %8.2f, V= %16.6f, enth= %f"%(Etol,P,V,enth)
            
        except:
            enth = 111111.11111111
            FinalPressure = 111111.11111111
            pass
        
        return enth,FinalPressure
