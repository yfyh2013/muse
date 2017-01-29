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

def Operator_Eval(Num_OperCalled,Num_OperSuccess,prob):
    FinalProbFile = open('../.FinalProb','w')
    
    for i in range(len(prob)):
        if Num_OperCalled[i] != 0:
            if float(Num_OperSuccess[i]) / float(Num_OperCalled[i]) >= 0.75:
                prob[i] += 0.1
            else:
                if prob[i] - 0.1 < 0.05:
                    prob[i] = 0.05
                else:
                    prob[i] -= 0.1
                    
    finalprob = [round(i / sum(prob),3) for i in prob]
    
    pickle.dump(finalprob,FinalProbFile)
    FinalProbFile.close()
    
    return finalprob
