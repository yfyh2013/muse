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

def Decompose2(N):
    i = 1
    allpair = []
    while i < N:
        allpair.append([i,N-i])
        i +=1
        
    allpair2 = []
    for a in allpair:
        if sorted(a) not in allpair2:
            allpair2.append(sorted(a))
    #print allpair2
    return allpair2


def Decompose3(N):
    i = 2
    allpair = []
    while i < N:
        for jj in Decompose2(i):
            allpair.append([])
            for kk in jj:
                allpair[-1].append(kk)
            allpair[-1].append(N-i)
        i +=1
    
    allpair2 = []
    for a in allpair:
        if sorted(a) not in allpair2:
            allpair2.append(sorted(a))
    #print allpair2
    return allpair2    

def Decompose4(N):
    i = 2
    allpair = []
    while i < N:
        for p1 in Decompose2(i):
            for p2 in Decompose2(N-i):
                allpair.append(p1+p2)
        i +=1
        
    allpair2 = []
    for a in allpair:
        if sorted(a) not in allpair2 and a != []:
            allpair2.append(sorted(a))
            
#    print allpair2
    return allpair2

def Decompose5(N):
    i = 1
    allpair = []
    while i < N:
        for p1 in Decompose4(N-i):
            allpair.append(p1+[i])
        i +=1
        
    allpair2 = []
    for a in allpair:
        if sorted(a) not in allpair2 and a != []:
            allpair2.append(sorted(a))
            
    return allpair2
    
def FindRatio(a,b):
    '''
    Find atom number ratio in Wyckoff positions
    at given space group number.
    '''
#    print a,b
    ratiolist = []
    used = []
    itema = []
    aa = []
    occ = []
    for ia,i in enumerate(a):
#        print "="*30
#        print "ia0:",ia,i
        success = False
        if i in b:
            for j,ib in enumerate(b):
                if i == ib and j not in used and itema.count(i) < a.count(i):
                    aa.append(i)
                    occ.append({indict['NameOfAtoms'][ia]:i})
                    used.append(j)
                    itema.append(i)
                    success = True
                    break
                else:
                    success = False

        if i not in b or not success:
            for ik,ii in enumerate(Decompose2(i)):
#                print ii
                a2 = []
                used2 = []
                for jj in ii:
                    if jj in b:
                        for j,ib in enumerate(b):
                            if jj == ib and j not in used and j not in used2:
                                used2.append(j)
                                a2.append(jj)
                                break

                if len(used2) == 2:
                    success = True
                    #print success
                    for ij,kk in enumerate(used2):
                        used.append(kk)
                        aa.append(a2[ij])
                        occ.append({indict['NameOfAtoms'][ia]:a2[ij]})
                    break
            
            if not success:
                for ik,ii in enumerate(Decompose3(i)):
                    a3 = []
                    used3 = []
                    for jj in ii:
                        if jj in b:
                            for j,ib in enumerate(b):
                                if jj == ib and j not in used and j not in used3:
                                    used3.append(j)
                                    a3.append(jj)
                                    break
                        
                    if len(used3) == 3:
                        success = True
                        for ij,kk in enumerate(used3):
                            used.append(kk)
                            aa.append(a3[ij])
                            occ.append({indict['NameOfAtoms'][ia]:a3[ij]})
                        break

#        print occ
                        
#            if not success:
#                for ik,ii in enumerate(Decompose4(i)):
#                    a4 = []
#                    used4 = []
#                    for jj in ii:
#                        if jj in b:
#                            for j,ib in enumerate(b):
#                                if jj == ib and j not in used and j not in used4:
#                                    used4.append(j)
#                                    a4.append(jj)
#                                    break
#        
#                    if len(used4) == 4:
#                        success = True
#                        for ij,kk in enumerate(used4):
#                            used.append(kk)
#                            aa.append(a4[ij])
#                            occ.append({indict['NameOfAtoms'][ia]:a4[ij]})
#                        break                        
#
#            if not success:
#                for ik,ii in enumerate(Decompose5(i)):
#                    a5 = []
#                    used5 = []
#                    for jj in ii:
#                        if jj in b:
#                            for j,ib in enumerate(b):
#                                if jj == ib and j not in used and j not in used5:
#                                    used5.append(j)
#                                    a5.append(jj)
#                                    break
#        
#                    if len(used5) == 5:
#                        success = True
#                        for ij,kk in enumerate(used5):
#                            used.append(kk)
#                            aa.append(a5[ij])
#                            occ.append({indict['NameOfAtoms'][ia]:a5[ij]})
#                        break    
#                        
        for jj in used:
            if jj not in ratiolist:
                ratiolist.append(jj)
    
    
#    print occ,aa,ratiolist
#    print sum(a),sum(aa)
    issuccess = False
    if sum(a) == sum(aa):
        issuccess = True
#    print occ,aa,ratiolist,issuccess
    return occ,aa,ratiolist,issuccess
                
if __name__ == '__main__':
    a = [64,64,64]
    b = [16,32,32,48,64,96,96,192]
    occ,aa,ratiolist,issuccess = FindRatio(a,b)
#    print Decompose2(96)
    print "="*30
    print a
    print aa
    print b
    print occ
    print ratiolist
        
