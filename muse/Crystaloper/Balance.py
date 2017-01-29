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

def Balance(randstr0,nnn):
    cl = randstr0.get_cell()
    Len_abc=[]
    for i in range(3):
        Len_abc.append((cl[i][0]**2+cl[i][1]**2+cl[i][2]**2)**.5)

    maxi = Len_abc.index(max(Len_abc))
    mini = Len_abc.index(min(Len_abc))

    if Len_abc[0] != Len_abc[1] and Len_abc[1] != Len_abc[2] and Len_abc[2] != Len_abc[1]:
        n1 = 0; n2 = 0; n3 = 0
        if maxi == 0: n1 = nnn[0]
        elif maxi == 1: n2 = nnn[0]
        elif maxi == 2: n3 = nnn[0]

        if mini == 0: n1 = nnn[2]
        elif mini == 1: n2 = nnn[2]
        elif mini == 2: n3 = nnn[2]

        if n1 == 0: n1 = nnn[1]
        elif n2 == 0: n2 = nnn[1]
        elif n3 == 0: n3 = nnn[1]

    elif Len_abc[0] == Len_abc[1] and Len_abc[1] == Len_abc[2] and Len_abc[2] == Len_abc[1]:
        n1,n2,n3 = nnn[0],nnn[1],nnn[2]

    elif Len_abc[0] == Len_abc[1] or Len_abc[1] == Len_abc[2] or Len_abc[2] == Len_abc[1]:
        n1 = 0; n2 = 0; n3 = 0
        if maxi == 0: n1 = nnn[0]
        elif maxi == 1: n2 = nnn[0]
        elif maxi == 2: n3 = nnn[0]

        if mini == 0: n1 = nnn[2]
        elif mini == 1: n2 = nnn[2]
        elif mini == 2: n3 = nnn[2]

        if n1 == 0: n1 = nnn[1]
        elif n2 == 0: n2 = nnn[1]
        elif n3 == 0: n3 = nnn[1]

    return n1,n2,n3
