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

def MakeInputfile(infile,outfile,tag,n):
    INfile=open(infile,'r')
    OUTfile=open(outfile,'w')
    line = INfile.readline()
    
    i= 1
    while line:
        if len(line.strip().split()) > 1 and tag in line.strip().split():
            
            while line and i == n:
                OUTfile.write(line)
                line = INfile.readline()
                if len(line.strip().split()) > 1 and tag in line.strip().split():
                    break
            i +=1

            line = INfile.readline()
        
        line = INfile.readline()
        
    INfile.close()
    OUTfile.close()
    
if __name__ == '__main__':
    #MakeInputfile('sys.cell','s.cell','#castepinput',4)
    MakeInputfile('INCARs','INCAR','Opt-step',3)
