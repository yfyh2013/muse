#! /bin/bash

rm -rf strana sym 2>/dev/null
mkdir strana

if [ $# -lt 2 ]; then
    echo "USAGE: strana.sh poscars-n_file accuracy"
    echo " e.g.: strana.sh poscars-2 0.2"
    exit
    fi

poscarfile=$1
Tolerance=$2
bestenthalpyfile=best-enthalpy-3
usefindsym=0   # 1 for yes, 0 for no
usespglib=1   # 1 for yes, 0 for no

cat > ana.py <<EOF
#! /usr/bin/python

import os
tolerance = $Tolerance

bestpos_file = open("$poscarfile","r")
line = bestpos_file.readline()

sname=line.split()[0]
nStr = 0
nstrs = 0
while line:                           # Determine # of structures
    nStr = line.count("%s"%sname)
    line = bestpos_file.readline()
    nstrs += int(nStr)
    
bestpos_file.seek(0)
i=1
tag=[]
tag_file = open("strana/tagfile","w")
tag_file.write(str(nstrs)+" ")
while i <= nstrs:                     # Generate findsym.in files
    line = bestpos_file.readline()
    line1 = line.strip().split()
    strtag = line1[1]
    filefindsym = "strana/findsym"+line1[1]+".in"
    findsymin_file =open(filefindsym,"w")
    
    findsymin_file.write(line)
    tag_file.write(str(strtag)+" ")
    findsymin_file.write(str(tolerance)+"\n")
    findsymin_file.write("1"+"\n")
    line = bestpos_file.readline()
    j = 0
    while j <= 2:
        line = bestpos_file.readline().strip().split()
        findsymin_file.write("%17.13f %17.13f %17.13f" %(float(line[0]),		float(line[1]),float(line[2])) + "\n")
        j += 1
        natomtype = []
        
    j=0
    line = bestpos_file.readline().strip().split()
    line = bestpos_file.readline().strip().split()
    natoms = 0
    while j < len(line):
        natomtype.append(line[j])
        natoms += int(natomtype[j])
        j += 1

    findsymin_file.write("1"+"\n")
    findsymin_file.write("1 0 0"+"\n")
    findsymin_file.write("0 1 0"+"\n")
    findsymin_file.write("0 0 1"+"\n")
    findsymin_file.write(str(natoms)+"\n")
    
    j=0
    while j < len(line):
        l = 1
        while l <= int(natomtype[j]):
            findsymin_file.write(str(j+1)+" ")
            l += 1
        j += 1

    findsymin_file.write("\n")
    line = bestpos_file.readline()

    j=1
    while j <= natoms:
        line = bestpos_file.readline().strip().split()
        findsymin_file.write("%17.13f %17.13f %17.13f" %(float(line[0]),		float(line[1]),float(line[2])) + "\n")
        j += 1
        
    i += 1
    findsymin_file.close()

bestpos_file.seek(0)
i=1
while i <= nstrs:                     # Generate POSCAR files
    
    j = 1
    while j <= 8 + natoms:
        line = bestpos_file.readline()
        if j == 1:
            fileposcar = "strana/POSCAR"+line.split()[1]+".vasp"
            poscar_file = open(fileposcar,"w")
        poscar_file.write(line)
        j += 1
    poscar_file.close()
    i += 1
bestpos_file.close()

EOF

python ana.py
nstrs=`cat strana/tagfile | awk '{print $1}'`

echo "Please wait..."

i=2
while [ $i -le $(($nstrs+1)) ]; do

n=`cat strana/tagfile | awk -v nvar="$i"  '{print $nvar}'`

Hmin=`tail -1 $bestenthalpyfile | awk '{print $4}'`

H=`grep " $n " $poscarfile | awk '{print $5}'`
deltaH=`echo $H $Hmin | awk '{printf"%.6f",1000*($1-$2)}'`


if [ $usefindsym = 1 ]; then

sg1=`findsym <findsym$n.in | grep "Space Group" | awk '{print $3}'`
sg2=`findsym <findsym$n.in | grep "Space Group" | awk '{print $4}'`
sg3=`findsym <findsym$n.in | grep "Space Group" | awk '{print $5}'`

fi


if [ $usespglib = 1 ]; then

cat >spg.py <<EOF
#! /usr/local/bin/python2.7

import sys,Generate_Supercell,Read_Write
from pyspglib import spglib

POS = Read_Write.read_vasp("strana/POSCAR"+"$n"+".vasp")
spg = spglib.get_spacegroup(POS,symprec=$Tolerance).split()

print spg[0],spg[1][1:-1]
    
EOF

sg4=`python spg.py | awk '{print $2}'`
sg5=`python spg.py | awk '{print $1}'`

fi

printf "%6s %12.4f %12.2f %8s %12s %8s %12s %12s\n" $n: $H $deltaH $sg4 $sg5 $sg1 $sg2 $sg3 >>sym

i=$(($i+1))
done

sort -k 2 -g -r sym > strinfo

cat strinfo
