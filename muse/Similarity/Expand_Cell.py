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

import copy
import random
from muse.Readwrite.Read_Write import *
from muse.Similarity.Charact_Triangles import *


def write_dict(chem, i, j, dd, posj, distang):
    if dd not in distang[i].keys():
        distang[i][dd] = [chem[j], list(posj)]

    elif dd in distang[i].keys() and list(posj) not in distang[i].values():
        while dd in distang[i].keys():
            dd += 10**(-11)*random.randint(0, 100)
        distang[i][dd] = [chem[j], list(posj)]
    
    return distang


def expand_cell(pc, distacc, angacc, r_cutoff):
    cell = pc.get_cell()
    pos = pc.get_scaled_positions()
    chem = pc.get_chemical_symbols()

    distang = {}
    for i, posi in enumerate(pos):
        distang[i] = {}
        posi0 = copy.deepcopy(posi)
        distang[i][0] = [chem[i], list(posi)]
        # in primitive cell
        for j, posj in enumerate(pos):
            if i != j:
                dd = dist_in_crycoord(copy.deepcopy(cell), posi, posj)
                distang = write_dict(chem, i, j, dd, posj, distang)
                # print j, dd, posi, posj, chem[j]
                # print distang
                # exit(0)

        # out primitive cell
        for j, posj in enumerate(pos):
            # shift along single direction
            posj0 = copy.deepcopy(posj)
            for k in range(3):
                for shift in (p1, m1):
                    posj = shift(posj0, k)
                    dd = dist_in_crycoord(copy.deepcopy(cell), posi0, posj)
                    distang = write_dict(chem, i, j, dd, posj, distang)

            # shift along two directions
            shifttwo1 = [[0, 1], [0, 2], [1, 2]]
            shifttwo2 = [[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]]
            
            howtwo1 = [[p1, p1], [m1, m1]]
            howtwo2 = [[p1, m1]]
            for ff in howtwo1:
                for jj in shifttwo1:
                    posj = ff[0](posj0, jj[0])
                    posj = ff[1](posj, jj[1])
                    dd = dist_in_crycoord(copy.deepcopy(cell), posi0, posj)
                    distang = write_dict(chem, i, j, dd, posj, distang)
                    
            for ff in howtwo2:
                for jj in shifttwo2:
                    posj = ff[0](posj0, jj[0])
                    posj = ff[1](posj, jj[1])
                    dd = dist_in_crycoord(copy.deepcopy(cell), posi0, posj)
                    distang = write_dict(chem, i, j, dd, posj, distang)
                        
            # shift along 3 directions
            howthree = [[p1, p1, p1], [p1, p1, m1], [p1, m1, m1], [m1, p1, p1],
                        [m1, m1, p1], [m1, p1, m1], [p1, m1, p1], [m1, m1, m1]]

            for ff in howthree:
                for jj in howthree:
                    posj = ff[0](posj0, 0)
                    posj = ff[1](posj, 1)
                    posj = ff[2](posj, 2)
                    dd = dist_in_crycoord(copy.deepcopy(cell), posi0, posj)
                    distang = write_dict(chem, i, j, dd, posj, distang)

    nearest = {}
    for i in range(len(distang.keys())):
        da = distang[i].keys()
        da.sort()
        nearest[i] = da

    cndistang = {}
    for natom in range(len(distang.keys())):
        # cndistang["%d-%s"%(natom+1,chem_sym[natom])] = {}
        # nearestdist = nearest[natom][1]
        ndistlist = []
        for dist in nearest[natom]:
            if abs(float(dist)) <= r_cutoff:
                ndistlist.append(float("%.16f" % dist))
        
        # print "Number of nearest neighbours is: %d"%len(ndistlist)
        # cndistang["%d-%s"%(natom+1,chem_sym[natom])]['cn'] = len(ndistlist)
        # cndistang["%d-%s"%(natom+1,chem_sym[natom])]['dist'] = ndistlist[0]
        # print natom, ndistlist

        if len(ndistlist) == 1:
            # print "Only one nearest neighbour!"
            allangleslist = []
        else:
            allangleslist = charact_triangles(cell, ndistlist, distang, natom, distacc, angacc)

        if allangleslist:
            cndistang["%d-%s" % (natom+1, chem[natom])] = allangleslist

    return cndistang

if __name__ == '__main__':
    pc1 = read_vasp('poscars-2-13-03.vasp')
    cnd1 = Nearest_Triangles(pc1)
    print cnd1
    print
    pc2 = read_vasp('poscars-2-12-16.vasp')
    cnd2 = Nearest_Triangles(pc2)
    print cnd2
