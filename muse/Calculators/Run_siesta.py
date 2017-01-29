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
from muse.Genstrs import Generate_Supercell
from muse.Readwrite import Read_Write,MakeInputfile
from muse.Readwrite.ReadInput import indict

def Run_siesta(nu,ng,nn,ParSubCom):
    Nopt = int(os.popen("grep Opt-step ../../siesta.fdf | wc -l").readline())
    
    atomorder = {indict['NameOfAtoms'][i]:\
                 [nu*int(indict['Num_atomsInChemForm'][i]),\
                  int(indict['AtomicNumber'][i])] for i in range(len(indict['NameOfAtoms']))}.items()
                  
    atomorder.sort()
    os.system("cp ../../*psf .")
    
    atoms = Read_Write.read_vasp("POSCAR")
    Read_Write.write_siesta((indict['NameSys'][0])+".STRUCT_IN",atoms)
    
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
        
        with open('sys.fdf','w') as includefile:
            print >>includefile, "SystemName  %12s"%(indict['NameSys'][0])
            print >>includefile, "SystemLabel  %11s\n"%(indict['NameSys'][0])
            print >>includefile, "NumberOfSpecies  %7d"%(len(indict['Num_atomsInChemForm']))
            print >>includefile, "NumberOfAtoms  %9d\n"%(nu*sum([int(x) for x in indict['Num_atomsInChemForm']]))
            
            print >>includefile, "%block ChemicalSpeciesLabel"
            ii=1; j=0
            for num in atoms.get_atomic_numbers():
                if j == 0: print >>includefile,"%2d %4d %6s"%(ii,num,Generate_Supercell.ChemicalSymbols[num])
                if j > 0:
                    if num != atoms.get_atomic_numbers()[j-1]:
                        ii +=1
                        print >>includefile,"%2d %4d %6s"%(ii,num,Generate_Supercell.ChemicalSymbols[num])
                        
                j +=1
            print >>includefile, "%endblock ChemicalSpeciesLabel\n"
            print >>includefile, "%block kgrid_Monkhorst_Pack"
            print >>includefile, "%2d %2d %2d %2d"%(K[0],0,0,0)
            print >>includefile, "%2d %2d %2d %2d"%(0,K[1],0,0)
            print >>includefile, "%2d %2d %2d %2d"%(0,0,K[2],0)
            print >>includefile, "%endblock kgrid_Monkhorst_Pack"
            
        # Write files and running.
        MakeInputfile.MakeInputfile('../siesta.fdf','siesta.fdf','Opt-step',i)
        
        os.system("sed -i 's/From-muse.in/%.2f/g' siesta.fdf"%(float(indict['Pressure'][0])))
        
        start = datetime.datetime.now()
        go = subprocess.Popen(ParSubCom,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        while go.poll() is None:
            time.sleep(2)
            now = datetime.datetime.now()
            if (now - start).seconds > float(indict['MaxHour'][0])*3600:
                os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
                break
            
        os.system("cp %s.STRUCT_OUT %s.STRUCT_IN"%(indict['NameSys'][0],indict['NameSys'][0]))
        os.system("cp OUT OUT%s"%(str(i)))
        os.system("rm INPUT_TMP*")
        with open('../log.muse','a') as logfile: print >>logfile, "Opt. %d ... done."%int(i)
        
        i += 1
        
    New_cry = Read_Write.read_siesta("%s.STRUCT_OUT"%(indict['NameSys'][0]),nu)
    Read_Write.write_vasp("CONTCAR",New_cry,label=indict['NameSys'][0]+": "+str(ng)+'-'+str(nn),\
                          direct=True,sort=True,vasp5=True)
                          
    try:
        enth = float(os.popen("grep eV/cell OUT | grep Free | tail -1 | awk '{print $5}'").readline())
        enth = float("%.8f"%enth)
        FinalPressure = float(os.popen("grep Pressure OUT | grep total | awk '{print $4/10}'").readline())
        #os.system("rm OUT")
        
    except:
        enth = 111111.11111111
        FinalPressure = 111111.11111111
        pass
    
    return enth,FinalPressure
