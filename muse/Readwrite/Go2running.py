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

import os
from muse.Readwrite.ReadInput import indict

muse = \
'''
+---------------------------------------------------------+
|             ###  #### ##    ## ################         |
|            #### #### ##    ## ##     #########          |
|           ##  ## ## ##    ## ##     ##                  |
|          ##   # ## ##    ## ###### #######              |
|         ##     ## ##    ##     ## ##                    |
|        ##     ## ########     ## ##                     |
|       ##     ################## #########               |
|                                                         |
|                      Ver. 2.2.0                         |
|                                                         |
|                            Let the atoms fly...         |
+---------------------------------------------------------+
'''

def Go2running():
    if not int(indict['IfPickUp'][0]):
        if os.path.isdir("running"):
            r=1
            while r < 1000:
                if not os.path.isdir("done%d"%r):
                    os.system("mv running done%d"%r)
                    break
                r +=1

        os.mkdir("running")

    if not os.path.isdir("running"): os.mkdir("running")

    os.chdir("running")
    os.system("cp ../muse.in .")

    os.mkdir("optdir")
    os.chdir("optdir")
    with open('../log.muse','w') as logfile:
        print >>logfile, muse

    return {},[]
