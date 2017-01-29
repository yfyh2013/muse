#! /usr/bin/python

from Readwrite import Read_Write
import os,sys
from pyspglib import spglib


if len(sys.argv) < 2:
    print
    print "USAGE: reana.py poscars-n_file accurcy"
    print " e.g.: reana.py poscars-4 0.2"
    print
    exit(0)

elif len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
    acc = 0.2
    print
    print "Using default accuracy: 0.2"
    print
    
elif len(sys.argv) == 3:
    acc = float(sys.argv[2])

os.system("mk_single_poscar.py %s"%sys.argv[1])

pos_file = open(sys.argv[1],"r")
line = pos_file.readline()
sname = line.split()[0]

enth = []
line = 1
while line:
    line = pos_file.readline()
    if line.count("%s"%sname):
        enth.append(float(line.split()[4]))


bigdict = {}
pos_file.seek(0)
line = 1
while line:
    line = pos_file.readline()
    if line.count("%s"%sname):
        fileposcar = "reana/POSCAR"+line.split()[1]+".vasp"
        POS = Read_Write.read_vasp(fileposcar)
        spg = spglib.get_spacegroup(POS,symprec=acc).split()

        bigdict[(float(line.split()[4])-min(enth))*1000] = [line.split()[1], spg[0], spg[1][1:-1], float(line.split()[4])] 

logfile = open("ordered-H"+"-"+sys.argv[1]+"-"+str(acc),'w')
sortedkeys = bigdict.keys()

sortedkeys.sort(reverse=True)
for i in sortedkeys:
    print >>logfile, "%5s: %14s %8s %16.6f %16.6f"%(bigdict[i][0],bigdict[i][1],bigdict[i][2],bigdict[i][3],i)
    print "%5s: %14s %8s %16.6f %16.6f"%(bigdict[i][0],bigdict[i][1],bigdict[i][2],bigdict[i][3],i)

print >>logfile, "="*70
print >>logfile, "ng-nn:             HM       sg            H(eV)      deltaH(meV)"
print >>logfile, "accuarcy = %s"%str(acc)
print >>logfile, "Deleted some strs. with the same enthalpy!!!"

print "="*70
print "ng-nn:             HM       sg            H(eV)      deltaH(meV)"
print "accuarcy = %s"%str(acc)
print "Deleted some strs. with the same enthalpy!!!"

print "Results saved to ordered-H"+"-"+sys.argv[1]+"-"+str(acc)
print

