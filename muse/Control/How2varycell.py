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
from muse.Operators.Crystal import Mut_Gauss,Random_Move,CrossOver,ExchAtoms,Ripple,Slip,Twist

def How2varycell(picked_operator,nu,ng,nn,Old_cry,BigDict,picked_indi,picked_indj,V_best,Num_OperSuccess,kT):
    if picked_operator == 'Mut':
        Op = Mut_Gauss.MUTGAUSS(nu,ng,nn,copy.deepcopy(Old_cry),\
                                BigDict,picked_indi,V_best,\
                                Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Mut()

    elif picked_operator == 'Cross':
        Op = CrossOver.CROSSOVER(nu,ng,nn,copy.deepcopy(Old_cry),\
                                 BigDict,picked_indi,picked_indj,\
                                 V_best,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Cross()

    elif picked_operator == 'CrossMut':
        Op = CrossOver.CROSSOVER(nu,ng,nn,copy.deepcopy(Old_cry),\
                                 BigDict,picked_indi,picked_indj,\
                                 V_best,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_CrossMut()

    elif picked_operator == 'Slip':
        Op = Slip.SLIP(nu,ng,nn,copy.deepcopy(Old_cry),BigDict,\
                       picked_indi,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Slip()

    elif picked_operator == 'Exch':
        Op = ExchAtoms.EXCHATOMS(nu,ng,nn,copy.deepcopy(Old_cry),\
                                 BigDict,V_best,picked_indi,\
                                 picked_indj,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Exch()

    elif picked_operator == 'ExchMut':
        Op = ExchAtoms.EXCHATOMS(nu,ng,nn,copy.deepcopy(Old_cry),\
                                 BigDict,V_best,picked_indi,\
                                 picked_indj,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_ExchMut()

    elif picked_operator == 'CrossExch':
        Op = ExchAtoms.EXCHATOMS(nu,ng,nn,copy.deepcopy(Old_cry),\
                                 BigDict,V_best,picked_indi,\
                                 picked_indj,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_CrossExch()

    elif picked_operator == 'RandomMove':
        Op = Random_Move.RANDOMMOVE(nu,ng,nn,copy.deepcopy(Old_cry),\
                                    BigDict,picked_indi,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Random_Move()

    elif picked_operator == 'Ripple':
        Op = Ripple.RIPPLE(nu,ng,nn,copy.deepcopy(Old_cry),BigDict,\
                           V_best,picked_indi,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Ripple()

    elif picked_operator == 'Twist':
        Op = Twist.TWIST(nu,ng,nn,copy.deepcopy(Old_cry),BigDict,\
                         picked_indi,Num_OperSuccess,kT)
        nk,Num_OperSuccess,BigDict = Op.Run_Twist()
        
    return nk,Num_OperSuccess,BigDict
