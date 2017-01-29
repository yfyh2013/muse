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

import numpy as np
import math


def calc_angle(l1, l2):
    v10 = l1
    v12 = l2
    v10 /= np.linalg.norm(v10)
    v12 /= np.linalg.norm(v12)
    angle = np.vdot(v10, v12)
    angle = np.arccos(angle)*180/math.pi

    return angle


def calc_3atoms_angle(pos1, pos2, pos3):
    v10 = pos2 - pos1
    v12 = pos3 - pos1
    v10 /= np.linalg.norm(v10)
    v12 /= np.linalg.norm(v12)
    angle = np.vdot(v10, v12)
    angle = np.arccos(angle)*180/math.pi

    return angle


def dihedral_angle(cell_a, cell_b, cell_c):
    a = cell_a
    b = cell_b
    c = cell_c

    axb = np.cross(a, b)
    axb /= np.linalg.norm(axb)

    bxc = np.cross(b, c)
    bxc /= np.linalg.norm(bxc)

    angle = np.vdot(axb, bxc)

    dangle = 180-np.arccos(angle)*180/math.pi
    # print dangle

    return dangle

if __name__ == '__main__':
    cell = np.array([[2.7085009435849550,   -2.7085009435849550,   -0.0000000000000000],
                     [-2.7085009435849550,    0.0000000000000000,   -2.7085009435849550],
                     [2.7085009435849550,    2.7085009435849550,   -0.0000000000000000]])

    cell_a = cell[0]
    cell_b = cell[1]
    cell_c = cell[2]

    dihedral_angle(cell_a, cell_b, cell_c)
    dihedral_angle(cell_b, cell_c, cell_a)
    dihedral_angle(cell_c, cell_a, cell_b)

