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

import pickle
from muse.Readwrite.ReadInput import indict

def Count_SelectTime(picked_operator,Num_OperCalled,ng):
    if len(indict['AtomicNumber']) > 1:
        if indict['Algorithm'][0] == 'MAC':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1
            elif picked_operator == 'CrossMut': Num_OperCalled[ng][2] +=1
            elif picked_operator == 'Slip': Num_OperCalled[ng][3] +=1
            elif picked_operator == 'RandomMove': Num_OperCalled[ng][4] +=1
            elif picked_operator == 'Ripple': Num_OperCalled[ng][5] +=1
            elif picked_operator == 'Twist': Num_OperCalled[ng][6] +=1
            elif picked_operator == 'Exch': Num_OperCalled[ng][7] +=1
            elif picked_operator == 'ExchMut': Num_OperCalled[ng][8] +=1
            elif picked_operator == 'CrossExch': Num_OperCalled[ng][9] +=1
            
        elif indict['Algorithm'][0] == 'EA':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1

        elif indict['Algorithm'][0] == 'EA2':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1
            elif picked_operator == 'Slip': Num_OperCalled[ng][2] +=1
            elif picked_operator == 'Twist': Num_OperCalled[ng][3] +=1
            
        elif indict['Algorithm'][0] == 'SA':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'RandomMove': Num_OperCalled[ng][1] +=1
        
        elif indict['Algorithm'][0] == 'BH':
            if picked_operator == 'RandomMove': Num_OperCalled[ng][0] +=1

        
    else:
        if indict['Algorithm'][0] == 'MAC':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1
            elif picked_operator == 'CrossMut': Num_OperCalled[ng][2] +=1
            elif picked_operator == 'Slip': Num_OperCalled[ng][3] +=1
            elif picked_operator == 'RandomMove': Num_OperCalled[ng][4] +=1
            elif picked_operator == 'Ripple': Num_OperCalled[ng][5] +=1
            elif picked_operator == 'Twist': Num_OperCalled[ng][6] +=1
            
        elif indict['Algorithm'][0] == 'EA':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1
            
        elif indict['Algorithm'][0] == 'EA2':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'Cross': Num_OperCalled[ng][1] +=1
            elif picked_operator == 'Slip': Num_OperCalled[ng][2] +=1
            elif picked_operator == 'Twist': Num_OperCalled[ng][3] +=1

        elif indict['Algorithm'][0] == 'SA':
            if picked_operator == 'Mut': Num_OperCalled[ng][0] +=1
            elif picked_operator == 'RandomMove': Num_OperCalled[ng][1] +=1
        
        elif indict['Algorithm'][0] == 'BH':
            if picked_operator == 'RandomMove': Num_OperCalled[ng][0] +=1
        
    with open('.Num_OperCalled','w') as Num_OperCalledFile:
        pickle.dump(Num_OperCalled,Num_OperCalledFile)
        
    return Num_OperCalled
