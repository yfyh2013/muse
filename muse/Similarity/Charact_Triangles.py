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

from muse.Similarity.Calc_Angle import calc_angle
from muse.Similarity.Similarity import if_eq_triang
from math import cos, acos, degrees, radians
import copy


def dist_in_crycoord(cl, pos1, pos2):
    if list(pos1) != list(pos2):
        ax = cl[0][0]
        ay = cl[0][1]
        az = cl[0][2]
        bx = cl[1][0]
        by = cl[1][1]
        bz = cl[1][2]
        cx = cl[2][0]
        cy = cl[2][1]
        cz = cl[2][2]
        
        a = (ax**2+ay**2+az**2)**0.5
        b = (bx**2+by**2+bz**2)**0.5
        c = (cx**2+cy**2+cz**2)**0.5
        alpha = radians(calc_angle(cl[1],cl[2]))
        beta = radians(calc_angle(cl[2],cl[0]))
        gamma = radians(calc_angle(cl[0],cl[1]))
        d1 = a**2*(pos1[0]-pos2[0])**2+b**2*(pos1[1]-pos2[1])**2+c**2*(pos1[2]-pos2[2])**2
        d2 = 2*b*c*cos(alpha)*(pos1[1]-pos2[1])*(pos1[2]-pos2[2])
        d3 = 2*a*c*cos(beta)*(pos1[0]-pos2[0])*(pos1[2]-pos2[2])
        d4 = 2*a*b*cos(gamma)*(pos1[0]-pos2[0])*(pos1[1]-pos2[1])
        d = (d1 + d2 + d3 + d4)**0.5
    
    else:
        d = 0
    
    return d


def p1(p, i):
    pp = copy.deepcopy(p)
    pp[i] += 1
    return pp


def m1(p, i):
    pp = copy.deepcopy(p)
    pp[i] -= 1
    return pp


def charact_triangles(cell, ndistlist, distang, natom, distacc, angacc):
    triangles = []
    for i, d1 in enumerate(ndistlist):
        for j, d2 in enumerate(ndistlist):
            if list(distang[natom][d1][1]) != list(distang[natom][d2][1]):
                # print
                # print natom, distang[natom][0], distang[natom][d1], distang[natom][d2]
                # print
                if i < j:
                    a = dist_in_crycoord(copy.deepcopy(cell), distang[natom][0][1], distang[natom][d1][1])
                    b = dist_in_crycoord(copy.deepcopy(cell), distang[natom][0][1], distang[natom][d2][1])
                    c = dist_in_crycoord(copy.deepcopy(cell), distang[natom][d1][1], distang[natom][d2][1])
                    # print a, b, c, chem[natom]

                    if a != 0 and b != 0 and c !=0:
                        abc = [(b**2+c**2-a**2)/2/b/c, (a**2+c**2-b**2)/2/a/c, (b**2+a**2-c**2)/2/b/a]

                        for ii, z in enumerate(abc):
                            if abs(abs(z) - 1.0) <= 0.0001:
                                if z > 0:
                                    abc[ii] = 1.0
                                if z < 0:
                                    abc[ii] = -1.0

                        newlist = [(a, acos(abc[0]), distang[natom][d2][0]),
                                   (b, acos(abc[1]), distang[natom][d1][0]),
                                   (c, acos(abc[2]), distang[natom][0][0])]
                        newlist.sort()

                        ifexist = False
                        if triangles == []:
                            triangles.append(newlist)
                        else:
                            for it in triangles:
                                ifexist = if_eq_triang(newlist, it, distacc, angacc)
                                if ifexist:
                                    break

                            if not ifexist:
                                triangles.append(newlist)

    return triangles
