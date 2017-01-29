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

from muse.Readwrite import Read_Write

def Centering(Cry):
    pos = Cry.get_scaled_positions()
    x_all = [pos[i][0] for i in range(len(pos))]
    y_all = [pos[i][1] for i in range(len(pos))]
    z_all = [pos[i][2] for i in range(len(pos))]
    deltax = 0.5 - 0.5*(min(x_all)+max(x_all))
    deltay = 0.5 - 0.5*(min(y_all)+max(y_all))
    deltaz = 0.5 - 0.5*(min(z_all)+max(z_all))
    
    for i in range(len(pos)):
        pos[i] = pos[i] + [deltax,deltay,deltaz]
        
    Cry.set_scaled_positions(pos)
    
    return Cry
    
if __name__ == '__main__':
    Cry= Read_Write.read_vasp('POSCAR')
    Centering(Cry)
    Read_Write.write_vasp("p.vasp",Cry,vasp5=True)
    print Cry
