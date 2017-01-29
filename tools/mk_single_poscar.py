#! /usr/bin/python

import os
import sys
from Readwrite import Read_Write
from pyspglib import spglib

os.system("rm -r reana 2>/dev/null")
os.mkdir("reana")

if len(sys.argv) < 2:
    print
    print "USAGE: mk_single_poscar.py poscars-n_file"
    print " e.g.: mk_single_poscar.py poscars-4"
    print
    exit(0)

pos_file = open(sys.argv[1],"r")
line = pos_file.readline()
sname = line.split()[0]

pos_file.seek(0)
while line:
    line = pos_file.readline()
    if line.count("%s"%sname):
        fileposcar = "reana/POSCAR"+line.split()[1]+".vasp"
        poscar_file = open(fileposcar,"w")
        poscar_file.write(line)
    else:
        poscar_file.write(line)
    
