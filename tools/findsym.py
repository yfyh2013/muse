#! /usr/local/bin/python2.7

import ase
import sys
from ase import io
from ase.lattice.spacegroup import crystal
from ase.lattice.spacegroup import findsym

a = crystal(ase.io.read(sys.argv[1],format="vasp"))

d=findsym.findsym(a,tol=5e-2)
print d['spacegroup']
