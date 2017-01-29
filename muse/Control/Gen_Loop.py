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

import os,copy,pickle,First_Gen,Selector,How2varycell
from muse.Calculators import DirectOpt
from muse.Readwrite import Print_Sep
from muse.Crystaloper import How2rep,DelRepAndCompl
from muse.Opereval import Operator_Eval,Count_SelectTime
from muse.Readwrite.ReadInput import indict
#from muse.Control import Basin_Hopping

def Gen_Loop(nu,nulist,BigDict):

    Num_OperCalled ={}; Num_OperSuccess ={}; prob = 0; ng = 1; nn0 = 0
    kT = float(indict['InitTemp_anneal'][0]) * 1.380658e-23 / 1.60217733e-19

    if int(indict['IfPickUp'][0]) == 1:
        RestartFile = open('../.Restart','r')
        BigDict = pickle.load(RestartFile)
        RestartFile.close()

        ng = int(os.popen("tail -1 ../all-enthalpy-%s | awk '{print $1}'"%str(nu)).readline().split('-')[0])
        nn0 = int(os.popen("tail -1 ../all-enthalpy-%s | awk '{print $1}'"%str(nu)).readline().split('-')[1][:-2])-1

        if ng >= 3:
            Num_OperSuccessFile = open('../.Num_OperSuccess','r')
            Num_OperCalledFile = open('../.Num_OperCalled','r')
            FinalProbFile = open('../.FinalProb','r')

            Num_OperSuccess = pickle.load(Num_OperSuccessFile)
            Num_OperCalled = pickle.load(Num_OperCalledFile)
            prob = pickle.load(FinalProbFile)

            Num_OperSuccessFile.close()
            Num_OperCalledFile.close()
            FinalProbFile.close()

        else:
            Num_OperSuccess[ng] = [1 for _ in range(len(Selector.opitems))]
            Num_OperCalled[ng] = [1 for _ in range(len(Selector.opitems))]
            prob = Selector.initprob

############# for test ##################
#        Num_OperSuccess[ng] = [1 for _ in range(len(Selector.opitems))]
#        Num_OperCalled[ng] = [1 for _ in range(len(Selector.opitems))]
#        prob = Selector.initprob
############# for test ##################

        if not ng in BigDict[nu]:
            ng -=1

        if nn0 >= int(indict['PopSize'][0])-1:
            ng +=1
            nn0 = 0

    while ng <= int(indict['MaxGen'][0]):
        if ng == 1:
            with open('../.sg','w') as sgfile:
                print >> sgfile, 0
        elif os.path.isfile('../.sg'):
            os.remove('../.sg')

        kT *= 0.90
        BigDict[nu][ng] = Print_Sep.Print_Sep(nu, ng)

        IfRep, Maxnu, RepTimes = How2rep.How2rep(nu,nulist)
        Num_OperCalled[ng] = [0 for _ in range(len(Selector.opitems))]
        Num_OperSuccess[ng] = [0 for _ in range(len(Selector.opitems))]

        if ng == 1:
            nk, BigDict = First_Gen.First_Gen(BigDict, nu, ng, IfRep, Maxnu, RepTimes)

        else:
            # Direct opt. for the best strs. in the last gen.
            Old_cry, V_best = Selector.GenOrder(BigDict, nu, ng-1)
            if nn0 == 0:
                BigDict = DirectOpt.DirectOpt(BigDict,copy.deepcopy(Old_cry), nu, ng)

            while nn0 < int(indict['PopSize'][0])-int(indict['Num_Keep'][0]):
                nn = nn0 + int(indict['Num_Keep'][0]) + 1

                if ng >= 3 and nn0 == 0:
                    if ng == 3:
                        prob=Selector.initprob
                    prob = Operator_Eval.Operator_Eval(Num_OperCalled[ng-1], Num_OperSuccess[ng-1], prob)

                picked_operator = Selector.OperSelector(prob, nu, ng)
                picked_indi, picked_indj = Selector.IndiSelector(BigDict,nu,ng-1)
                Num_OperCalled = Count_SelectTime.Count_SelectTime(picked_operator,Num_OperCalled,ng)

                # picked_operator = 'CrossMut'  ## for test !!!!!
                nk, Num_OperSuccess,BigDict = How2varycell.How2varycell(picked_operator, nu, ng, nn, \
                                                          Old_cry, BigDict, picked_indi,\
                                                          picked_indj, V_best, Num_OperSuccess, kT)
                nn0 +=1

        if ng > 1 and int(indict['IfCompeteBetweenGens'][0]):
            Old_cry = DelRepAndCompl.Complement(BigDict, nu, ng, nn)

        # print ng
        if nn0 >= int(indict['PopSize'][0])-int(indict['Num_Keep'][0]):
            nn0 = 0
        ng = nk; ng += 1; #ns = ng + 1

    # Pure basin hopping, optionally
    # BigDict = Basin_Hopping.Basin_Hopping(BigDict,nu,ns,kT)

    with open('../log.muse', 'a') as logfile:
        print >>logfile
        print >>logfile, "All generations done."
        print >>logfile, "Good luck!"
        print >>logfile

    return BigDict
