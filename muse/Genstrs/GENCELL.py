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

class GENCELL:
    def __init__(self,sg_num):
        self.sg_num = sg_num
        self.cellp = []
        for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
        for i in range(3): self.cellp.append(45 + 90.0 * random.random())
        
    def max2min(self,cellp):
        listrue = 1.0 < max(cellp[0],cellp[1],cellp[2])/min(cellp[0],cellp[1],cellp[2]) < 5.0
        return listrue
    
    def isang(self,cellp):
        alpha,beta,gamma = cellp[-3:]
        alpha *= math.pi/180.0
        beta *= math.pi/180.0
        gamma *= math.pi/180.0
        cx = math.cos(beta)
        cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)
        isangtrue = cx*cx + cy*cy <= 1
        return isangtrue
        
    def ifabeq(self,a,b):
        while min(a,b) < 1.0:
            a *= 10
            b *= 10
        
        return abs(a-b) < 1.0
        
    def ifeq90(self,alpha):
        return abs(alpha-90) < 10
        
    def Gen_sg_latt(self):
        listrue = self.max2min(self.cellp)
        isangtrue = self.isang(self.cellp)
        
        # Triclinic 
        if 1 <= self.sg_num <= 2:
            listrue = 0
            while not listrue or not isangtrue >= 1:
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                for i in range(3): self.cellp.append(45 + 90 * random.random())
                isangtrue = self.isang(self.cellp)            
                listrue = self.max2min(self.cellp)
        
        # Monoclinic 
        elif 3 <= self.sg_num <= 15:
            listrue = 0
            while not listrue or self.ifabeq(self.cellp[0],self.cellp[1]) or \
                                 self.ifabeq(self.cellp[1],self.cellp[2]) or \
                                 self.ifabeq(self.cellp[2],self.cellp[0]) or \
                                 not isangtrue or self.ifeq90(self.cellp[5]):
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                self.cellp.append(90)
                self.cellp.append(90)
                self.cellp.append(45 +90 * random.random())
                listrue = self.max2min(self.cellp)
                isangtrue = self.isang(self.cellp)
        
        # Orthorhombic 
        elif 16 <= self.sg_num <= 74:
            listrue = 0
            while not listrue or self.ifabeq(self.cellp[0],self.cellp[1]) or \
                                 self.ifabeq(self.cellp[1],self.cellp[2]) or \
                                 self.ifabeq(self.cellp[2],self.cellp[0]):
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                self.cellp.append(90)
                self.cellp.append(90)
                self.cellp.append(90)
                listrue = self.max2min(self.cellp)
        
        # Tetragonal         
        elif 75 <= self.sg_num <= 142:
            listrue = 0
            while not listrue or self.ifabeq(self.cellp[1],self.cellp[2]):
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                self.cellp[0] = self.cellp[1]
                self.cellp.append(90)
                self.cellp.append(90)
                self.cellp.append(90)
                listrue = self.max2min(self.cellp)
        
        # Trigonal 
        elif 143 <= self.sg_num <= 167:
            self.cellp[1] = self.cellp[0]
            self.cellp[2] = self.cellp[0]
            self.cellp[4] = self.cellp[3]
            self.cellp[5] = self.cellp[3]
            isangtrue = self.isang(self.cellp)          
            while not isangtrue or self.ifeq90(self.cellp[3]):
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                self.cellp[1] = self.cellp[0]
                self.cellp[2] = self.cellp[0]
                for i in range(3): self.cellp.append(45 + 90 * random.random())
                self.cellp[4] = self.cellp[3]
                self.cellp[5] = self.cellp[3]
                isangtrue = self.isang(self.cellp)
            
        # Hexagonal
        elif 168 <= self.sg_num <= 194:
            listrue = 0
            while not listrue or self.ifabeq(self.cellp[1],self.cellp[2]):
                self.cellp = []
                for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
                self.cellp[0] = self.cellp[1]
                self.cellp.append(90)
                self.cellp.append(90)
                self.cellp.append(120)
                listrue = self.max2min(self.cellp)
        
        # Cubic 
        elif 195 <= self.sg_num <= 230:
            self.cellp = []
            for i in range(3): self.cellp.append(random.uniform(0.1,1.0))
            self.cellp[1] = self.cellp[0]
            self.cellp[2] = self.cellp[0]
            self.cellp.append(90)
            self.cellp.append(90)
            self.cellp.append(90)
    
        return self.cellp

if __name__ == '__main__':
    num = random.randint(1,230)
    c = GENCELL(num)
    
    print num,c.Gen_sg_latt()
