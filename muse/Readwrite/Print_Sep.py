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

def Print_Sep(nu,ng):
    with open('../log.muse','a') as logfile:
        if nu <10:
            print >>logfile, str(nu)*22,"("+"ng"+"%d"%(ng)+")",str(nu)*22
            print >>logfile, "="*51
        else:
            print >>logfile, str(nu)*11,"("+"ng"+"%d"%(ng)+")",str(nu)*11
            print >>logfile, "="*51
            
    return {}
