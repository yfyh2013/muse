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

from spglib import spglib

def Findspg(atoms,prec=None):
    if prec == None:
        from muse.Readwrite.ReadInput import indict
        prec = float(indict['SymPrec'][0])

    spg0 = spglib.get_spacegroup(atoms,symprec=prec)
    if spg0:
        spg1 = spg0.split()
        spg = [str(spg1[0]), int(spg1[1][1:-1])]
    else:
        spg = []
    # print spg0,spg

    return spg
    
if __name__ == '__main__':
    from muse.Readwrite import Read_Write
    pc = Read_Write.read_vasp('194.vasp')
    dataset = spglib.get_symmetry_dataset(pc, symprec=0.1)
    spg = Findspg(pc,prec=0.1)
    print spg,dataset
