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

from muse.Similarity.Evaluate_Similarity import evaluate_similarity


def picker(Big_Dict, New_Big_Dict, done_list, enth_list, dd, natoms, E_Diff_max, r_cutoff, distacc, angacc):
    for enth1 in dd:
        for enth2 in dd:
            if enth1 in Big_Dict and enth2 in Big_Dict and abs(enth1 - dd[0])/natoms <= E_Diff_max and \
                                    abs(enth2 - dd[0])/natoms <= E_Diff_max:
                if dd.index(enth1) < dd.index(enth2) and (enth1, enth2) not in done_list:
                    # print (enth1, enth2), done_list
                    done_list.append((enth1, enth2))

                    no1 = Big_Dict[enth1][0]
                    no2 = Big_Dict[enth2][0]
                    spg1 = Big_Dict[enth1][1]
                    spg2 = Big_Dict[enth2][1]
                    comres = evaluate_similarity(Big_Dict[enth1][2], Big_Dict[enth2][2], distacc, angacc, r_cutoff)

                    if comres and spg1 == spg2:
                        ifsame = "YES, with same symmetry."

                    elif comres and spg1 != spg2:
                        ifsame = "YES, without same symmetry."

                    elif not comres and spg1 == spg2:
                        ifsame = "NO, with same symmetry."
                        if enth2 == enth_list[-1] and enth2 in Big_Dict:
                            New_Big_Dict[enth2] = Big_Dict[enth2]

                    elif not comres and spg1 != spg2:
                        ifsame = "NO, without same symmetry."
                        if enth2 == enth_list[-1] and enth2 in Big_Dict:
                            New_Big_Dict[enth2] = Big_Dict[enth2]

                    print "Are", no1, "(", spg1, enth1, ")", "and", no2, "(", spg2, enth2, ") identical?", ifsame

                    if enth1 not in New_Big_Dict and enth1 in Big_Dict:
                        New_Big_Dict[enth1] = Big_Dict[enth1]

                    if comres and enth2 in Big_Dict:
                        del Big_Dict[enth2]

    return None
