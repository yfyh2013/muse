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

try:
    from setuptools import setup,find_packages
except ImportError:
    from distutils.core import setup

setup(
      name = "muse",
      version = "2.2.0",
      description = "Package for predicting crystal structures",
      author = "Zhong-Li Liu",
      url = "http://sourceforge.net/projects/pymuse/",
      license = "LGPL",
      package_dir = {'muse':'muse'},
      package_data = {'muse.Symmetry':['data/*.dat']},
      packages = ['muse',
                'muse.Calculators',
                'muse.Control',
                'muse.Crystaloper',
                'muse.Genstrs',
                'muse.Operators',
                'muse.Operators.Crystal',
                'muse.Opereval',
                'muse.Readwrite',
                'muse.Similarity',
                'muse.Symmetry'],
      scripts = ['muse/muse','tools/hor'],
      )
