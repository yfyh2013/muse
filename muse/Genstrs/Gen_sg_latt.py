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

import random,math

def Gen_sg_latt(sg_num):
    cellp = []

    for i in range(3): cellp.append(random.uniform(0.1,1.0))
    for i in range(3): cellp.append(45 + 90.0 * random.random())

    listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0
    alpha,beta,gamma = cellp[-3:]
    alpha *= math.pi/180.0
    beta *= math.pi/180.0
    gamma *= math.pi/180.0
    cx = math.cos(beta)
    cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)     

    if 1 <= sg_num <= 2:
        listrue = 0
        while not listrue or cx*cx + cy*cy >= 1:
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            for i in range(3): cellp.append(45 + 90 * random.random())
            
            alpha,beta,gamma = cellp[-3:]
            alpha *= math.pi/180.0
            beta *= math.pi/180.0
            gamma *= math.pi/180.0
            cx = math.cos(beta)
            cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)             
            listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0
            #print "Here:%d  %6.2f"%(sg_num,1-cx*cx - cy*cy)

    elif 3 <= sg_num <= 15:
        listrue = 0
        while not listrue or cellp[5] == 90 or cx*cx + cy*cy >= 1:
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            cellp.append(90)
            cellp.append(90)
            cellp.append(45 +90 * random.random())
            listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0
            alpha,beta,gamma = cellp[-3:]
            alpha *= math.pi/180.0
            beta *= math.pi/180.0
            gamma *= math.pi/180.0
            cx = math.cos(beta)
            cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)             
            #print "Here:%d  %6.2f"%(sg_num,1-cx*cx - cy*cy)

    elif 16 <= sg_num <= 74:
        listrue = 0
        while not listrue or cellp[0] == cellp[1] or cellp[1] == cellp[2] or cellp[2] == cellp[0]:
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            cellp.append(90)
            cellp.append(90)
            cellp.append(90)
            listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0
            
    elif 75 <= sg_num <= 142:
        listrue = 0
        while not listrue or cellp[1] == cellp[2]:
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            cellp[0] = cellp[1]
            cellp.append(90)
            cellp.append(90)
            cellp.append(90)
            listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0

    elif 143 <= sg_num <= 167:
        cellp[1] = cellp[0]
        cellp[2] = cellp[0]
        cellp[4] = cellp[3]
        cellp[5] = cellp[3]
        alpha,beta,gamma = cellp[-3:]
        alpha *= math.pi/180.0
        beta *= math.pi/180.0
        gamma *= math.pi/180.0
        cx = math.cos(beta)
        cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)          
        #print "Here:%d  %6.2f"%(sg_num,1-cx*cx - cy*cy)
        while cx*cx + cy*cy >= 1 or cellp[3] == 90:
            #print "I am here"
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            cellp[1] = cellp[0]
            cellp[2] = cellp[0]
            for i in range(3): cellp.append(45 + 90 * random.random())
            cellp[4] = cellp[3]
            cellp[5] = cellp[3]
            alpha,beta,gamma = cellp[-3:]
            alpha *= math.pi/180.0
            beta *= math.pi/180.0
            gamma *= math.pi/180.0
            cx = math.cos(beta)
            cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)              
            #print cx,cy
            #print "Here:%d  %6.2f"%(sg_num,1-cx*cx - cy*cy)

    elif 168 <= sg_num <= 194:
        listrue = 0
        while not listrue or cellp[1] == cellp[2]:
            cellp = []
            for i in range(3): cellp.append(random.uniform(0.1,1.0))
            cellp[0] = cellp[1]
            cellp.append(90)
            cellp.append(90)
            cellp.append(120)
            listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0

    elif 195 <= sg_num <= 230:
        cellp = []
        for i in range(3): cellp.append(random.uniform(0.1,1.0))
        cellp[1] = cellp[0]
        cellp[2] = cellp[0]
        cellp.append(90)
        cellp.append(90)
        cellp.append(90)

    return cellp

if __name__ == '__main__':
    c = Gen_sg_latt(random.uniform(1,230))
    print c
