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

def Splitnum(nu):
    N = nu
    i = 2
    listnum = []
    while i < N:
        if N % i == 0:
            N /= i
            listnum.append(i)
            i -=1
        i +=1
    listnum.append(N)
    
    if len(listnum) > 3:
        while len(listnum) > 3:
            listnum.sort()
            x = listnum[0]*listnum[1]
            del(listnum[0])
            listnum[0] = x
            
    elif len(listnum) < 3:
        while len(listnum) != 3:
            listnum.append(1)
    listnum.sort()
    
    return listnum
    
if __name__ == '__main__':
    for _ in range(100):
        nnn = Splitnum(12)
        print nnn
