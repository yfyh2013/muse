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

from muse.Similarity import Build_Cell
from spglib import spglib


def reduce_cell(pc, distacc, angacc):
    primcell = spglib.find_primitive(pc, symprec=distacc, angle_tolerance=angacc)
    pcell = []
    for i in range(len(primcell[0])):
        pcell.append([primcell[0][i][j] for j in range(3)])

    pos = []
    for i in range(len(primcell[1])):
        pos.append([primcell[1][i][j] for j in range(3)])

    sys_name = [Build_Cell.ChemicalSymbols[ii] for ii in primcell[2]]
    cell_reduced = Build_Cell.build_cell(sys_name, cell=pcell, scaled_positions=pos)

    return cell_reduced
