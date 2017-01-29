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

import random


def decide_component(akeys):
    # judge component if the same
    name_dict = {}
    for namei in akeys:
        if namei.split("-")[1] not in name_dict:
            name_dict[namei.split("-")[1]] = 1
        else:
            name_dict[namei.split("-")[1]] += 1

    return name_dict


def if_eq_triang(distang1, distang2, distacc, angacc):
    iftrue = False

    # element name
    ifsameelement = distang1[0][2] == distang2[0][2] and \
                    distang1[1][2] == distang2[1][2] and \
                    distang1[2][2] == distang2[2][2]

    if not ifsameelement:
        iftrue = False

    else:
        # distance
        ifdisteq = abs(distang1[0][0]-distang2[0][0]) <= distacc and \
                   abs(distang1[1][0]-distang2[1][0]) <= distacc and \
                   abs(distang1[2][0]-distang2[2][0]) <= distacc

        if not ifdisteq:
            iftrue = False

        else:
            # angle
            anglist1 = [distang1[0][1], distang1[1][1], distang1[2][1]]
            anglist2 = [distang2[0][1], distang2[1][1], distang2[2][1]]
            ifangeq = abs(max(anglist1) - max(anglist2)) < angacc and \
                      abs(min(anglist1) - min(anglist2)) < angacc

            if ifangeq:
                iftrue = True

    return iftrue


def similarity(adict, bdict, distacc, angacc):
    # adict,bdict = bdict,adict
    a = adict.keys()
    b = bdict.keys()
    # print a
    # print b

    aname_dict = decide_component(a)
    bname_dict = decide_component(b)
    # print aname_dict, bname_dict

    ratios = []
    for name in aname_dict:
        if name in bname_dict:
            ratios.append(aname_dict[name] / bname_dict[name])

        else:
            print "Well, they have no same component; of course they are in different structures."
            # exit(0)

    # print len(aname_dict), len(bname_dict), ratios
    if ratios:
        if len(aname_dict) != len(bname_dict) or min(ratios) != max(ratios):
            print "Well, they have no same component; of course they are in different structures."
            # exit(0)

    # print aname_dict, bname_dict

    # calculate similarity of all the triangles formed with the reference atom
    similaritylist = []
    part_sim = []
    if len(a) == len(b):
        name_list = []
        aname_dict = {}
        bname_dict = {}

        for i, ai in enumerate(a):
            name = ai.split("-")[1]

            if name not in name_list:
                name_list.append(name)

        # print name_list

        for n, name in enumerate(name_list):
            aname_dict[name] = []
            for i, ai in enumerate(a):
                if ai.split("-")[1] == name:
                    aname_dict[name].append(ai)

        for n, name in enumerate(name_list):
            bname_dict[name] = []
            for j, bj in enumerate(b):
                if bj.split("-")[1] == name:
                    bname_dict[name].append(bj)
        len_list = []
        for it in name_list:
            aindex = range(len(aname_dict[it]))
            # print len(aname_dict[it])
            len_list.append([len(aname_dict[it]), it])

        len_list.sort()
        # print len_list
        for it in len_list:

            tot_sim = 0
            aindex = range(len(aname_dict[it[1]]))
            # bindex = list(itertools.permutations(range(len(bname_dict[it[1]]))))[:10]
            bindex = []

            i = 1
            num_try = 1
            while i <= len(bname_dict[it[1]]):
                num_try *= i
                i += 1
                if num_try > 100:
                    break

            while num_try:
                max_list = 1000
                while max_list:
                    if max_list == 1000:
                        ii = range(len(bname_dict[it[1]]))
                    else:
                        ii = random.sample(range(len(bname_dict[it[1]])), len(bname_dict[it[1]]))

                    if ii not in bindex:
                        bindex.append(ii)
                        break
                    max_list -= 1

                # print "ii=",ii
                sim = []
                for jj in range(len(aindex)):
                    # print aname_dict[it][aindex[jj]], bname_dict[it][list(ii)[jj]]
                    # print jj
                    aaa = adict[aname_dict[it[1]][aindex[jj]]]
                    bbb = bdict[bname_dict[it[1]][list(ii)[jj]]]
                    aaa.sort()
                    bbb.sort()
                    lena = len(aaa)
                    lenb = len(bbb)
                    # print aaa, bbb, lena, lenb
                    count_namei = 0
                    for i, ai in enumerate(aaa):
                        for j, aj in enumerate(bbb):
                            if bbb[j] is not None and if_eq_triang(aaa[i], bbb[j], distacc, angacc):
                                count_namei += 1
                                # print count_namei
                                # print count_namei
                                bbb[j] = None
                                break

                    simil = round(count_namei/float(min(lena, lenb)), 4)

                    sim.append(simil)

                # print sim
                # print sum(sim)
                if sum(sim) > tot_sim:
                    tot_sim = sum(sim)
                    part_sim = sim

                # print num_try
                if sum(sim) / float(len(sim)) == 1.0:
                    break

                num_try -= 1

            if part_sim != []:
                if sum(part_sim) / float(len(part_sim)) >= 0.9:
                    similaritylist.append(part_sim)
                else:
                    similaritylist.append([])
                    break
            # print similaritylist
        bindex = []
    # print bdict
    # print
    # print "similarity=", similaritylist
    # totalsimilarity = float(sum(similaritylist)) / float(len(similaritylist))
    # print totalsimilarity
    # print
    # print similaritylist
    # print sum(similaritylist)
    return similaritylist
