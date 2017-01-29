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
from muse.Similarity.Expand_Cell import *
from muse.Readwrite.Read_Write import *
from muse.Similarity.Similarity import similarity
from muse.Similarity.Reduce import reduce_cell
from spglib import spglib


def evaluate_similarity(pc1, pc2, distacc, angacc, r_cutoff):

    pc1_reduced = reduce_cell(pc1, distacc, angacc)
    pc2_reduced = reduce_cell(pc2, distacc, angacc)
    ifsame = False

    if int(pc1_reduced.get_number_of_atoms()) == int(pc2_reduced.get_number_of_atoms()):
        cnd1 = expand_cell(pc1_reduced, distacc, angacc, r_cutoff)
        cnd2 = expand_cell(pc2_reduced, distacc, angacc, r_cutoff)

        volratio = pc1.get_volume() / pc2.get_volume()
        if 0.90 <= volratio <= 1.10:
            pass
            # print "NOTE: Comparing at the same volume."
        else:
            print "NOTE: Comparing not at the same volume, the results may be further analysized."

        similaritylist = similarity(cnd1, cnd2, distacc, angacc)

        ncount = 0
        nn = 0
        for ss in similaritylist:
            for s in ss:
                nn += 1
                if s >= 0.9:
                    ncount += 1

        if similaritylist != [] and [] not in similaritylist and ncount >= 0.9 * nn:
            ifsame = True

    elif int(pc1_reduced.get_number_of_atoms()) != int(pc2_reduced.get_number_of_atoms()):
        # print "They are not identical due to different inequivalent atoms in the primitive cell."
        ifsame = False

    return ifsame

if __name__ == '__main__':
    pc1 = read_vasp('13-01.vasp')
    pc2 = read_vasp('06-02.vasp')
    Isthesame = CompDistTriang(pc1, pc2, 5)
    if Isthesame:
        print "They are the same structure."
    else:
        print "They are not the same structure."
