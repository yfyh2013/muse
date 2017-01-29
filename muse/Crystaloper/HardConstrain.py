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

from muse.Similarity import Calc_Angle
from muse.Readwrite.ReadInput import indict

def HardConstrain(O_cry,nu,Ang_min=35,Ang_max=145):
    total_atoms = O_cry.get_number_of_atoms()
    cl = O_cry.get_cell()
    ang1 = Calc_Angle.calc_angle(cl[0],cl[1])
    ang2 = Calc_Angle.calc_angle(cl[1],cl[2])
    ang3 = Calc_Angle.calc_angle(cl[2],cl[0])

    if min(ang1,ang2,ang3) >= Ang_min and max(ang1,ang2,ang3) <= Ang_max:
        L_angle=False
    else:
        L_angle=True

    cell = O_cry.get_cell()

    cell_a = cell[0]
    cell_b = cell[1]
    cell_c = cell[2]

    dangle1 = Calc_Angle.dihedral_angle(cell_a,cell_b,cell_c)
    dangle2 = Calc_Angle.dihedral_angle(cell_b,cell_c,cell_a)
    dangle3 = Calc_Angle.dihedral_angle(cell_c,cell_a,cell_b)

    #print dangle1,dangle2,dangle3
    if min(dangle1,dangle2,dangle3) >= 15 and max(dangle1,dangle2,dangle3) <= 165:
        L_dangle=False
    else:
        L_dangle=True

    Len_abc=[]
    for i in range(3):
        Len_abc.append((cell[i][0]**2+cell[i][1]**2+cell[i][2]**2)**.5)

    if min(Len_abc)>=float(indict['MinVect'][0]) and max(Len_abc) < 5*nu*min(Len_abc):
        L_len=False
    else:
        L_len=True

    adist=[]
    supercell222=O_cry.repeat((2,2,2))

    for i in range(total_atoms*8):
        for j in range(total_atoms*8):
            if j > i:
                adist.append(supercell222.get_distance(i,j))

    if min(adist)>=float(indict['MinDist'][0]):
        L_dist=False
    else:
        L_dist=True

    return L_dangle,L_angle,L_len,L_dist

