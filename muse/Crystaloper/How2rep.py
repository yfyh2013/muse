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

from muse.Readwrite.ReadInput import indict

def How2rep(nu,nulist):
    Maxnu = 0
    RepTimes = 0
    IfRep = False
    if nu > int(indict['MinNum_PrimitiveCell'][0]):
        IntTimes=[]
        for i in nulist:
            if max(nulist)%i == 0 and max(nulist) != i:
                IfRep = True
                IntTimes.append(i)
                
        if IntTimes != []:
            Maxnu = max(IntTimes)
            RepTimes = int(max(nulist)/max(IntTimes))
            
    return	IfRep,Maxnu,RepTimes
