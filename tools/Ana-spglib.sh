#! /bin/bash

rm -rf strana sym 2>/dev/null
mkdir strana

if [ $# -lt 2 ]; then
    echo "USAGE: Ana-spglib.sh poscars-n_file best-enthalpy-n_file accuracy"
    echo " e.g.: Ana-spglib.sh poscars-2 best-enthalpy-2 0.2"
    exit
    fi

poscarfile=$1
bestenthalpyfile=$2
Tolerance=$3

cat > mkposcar.py <<EOF
#! /usr/bin/python

import os
tolerance = $Tolerance

pos_file = open("$poscarfile","r")
line = pos_file.readline()
sname = line.split()[0]

pos_file.seek(0)
nstrs = 0
while line:                           # Determine # of structures
    line = pos_file.readline()
    nstrs += int(line.count("%s"%sname))
    
i=1
tag=[]
tag_file = open("strana/tagfile","w")
tag_file.write(str(nstrs)+" ")

pos_file.seek(0)

for i in range(7):
    line = pos_file.readline()

natoms=sum([int(item) for item in line.split()])

pos_file.seek(0)
i=1
while i <= nstrs:                   # Generate POSCAR files
    j = 1
    while j <= 8 + natoms:
        line = pos_file.readline()
        if j==1:
            tag_file.write(str(line.split()[1])+" ")
        if j == 1:
            fileposcar = "strana/POSCAR"+line.split()[1]+".vasp"
            poscar_file = open(fileposcar,"w")
        poscar_file.write(line)
        j += 1
    poscar_file.close()
    i += 1
pos_file.close()

EOF

python mkposcar.py
nstrs=`cat strana/tagfile | awk '{print $1}'`

echo
echo $nstrs structures in total.
echo
echo "Please wait..."

i=2
while [ $i -le $(($nstrs+1)) ]; do

n=`cat strana/tagfile | awk -v nvar="$i"  '{print $nvar}'`

Hmin=`tail -1 $bestenthalpyfile | awk '{print $4}'`

H=`grep " $n " $poscarfile | awk '{print $5}'`
deltaH=`echo $H $Hmin | awk '{printf"%.6f",1000*($1-$2)}'`

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

printf "%6s %12.4f %12.2f %8s %12s %8s %12s %12s\n" $n: $H $deltaH $sg4 $sg5 $sg1 $sg2 $sg3 >>sym

i=$(($i+1))
done

sort -k 2 -g -r sym > ordered-H-$3

cat ordered-H-$3

#rm mkposcar.py spg.py strana/tagfile
