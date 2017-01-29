#! /usr/bin/python

import os,sys
from Genstrs import Generate_Supercell
from Readwrite.Read_Write import read_vasp
from Readwrite.Read_Write import write_vasp
from pyspglib import spglib

if len(sys.argv) < 2:
    print
    print "USAGE: spg.py single-poscar_file accurcy"
    print " e.g.: spg.py 160.vasp 0.2"
    print
    exit(0)

elif len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
    acc = 0.2
    print
    print "Using default accuracy: 0.2"
    print
    
elif len(sys.argv) == 3:
    acc = float(sys.argv[2])

POS = read_vasp(sys.argv[1])
spg = spglib.get_spacegroup(POS,symprec=acc).split()
pc = spglib.find_primitive(POS,symprec=acc, angle_tolerance=-1.0)

print spg

if None not in pc:
    pcell = []
    for i in range(len(pc[0])):
        pcell.append([pc[0][i][j] for j in range(3)])
    
    pos = []
    for i in range(len(pc[1])):
        pos.append([pc[1][i][j] for j in range(3)])
    
    Ele = []
    for i in range(len(pc[2])):
        Ele.append(Generate_Supercell.ChemicalSymbols[pc[2][i]])
    
    primcell = Generate_Supercell.Generate_Supercell(Ele,cell=pcell,scaled_positions=pos)
    write_vasp('prim-'+sys.argv[1],primcell,label='prim-cell',direct=True,sort=True,vasp5=True)
else:
    print "This is already the primitive cell!!"
    
