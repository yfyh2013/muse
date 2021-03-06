#! /usr/local/bin/python2.7

import re,os,sys,math
import numpy as np
from numpy import matrix
from numpy import array
from ase.io.vasp import read_vasp
from pyspglib import spglib
from ase.lattice.spacegroup import cell

if len(sys.argv) < 3:
    print
    print "USAGE: pick_up_poscar.py poscars-x xx-xx"
    print " e.g.: pick_up_poscar.py poscars-2 03-21"
    print
    exit(0)

elif len(sys.argv) == 3:
    acc = 0.1
    print

elif len(sys.argv) == 4:
    acc = float(sys.argv[3])
    print

filename = sys.argv[1]
num = sys.argv[2]

pos_file = open("%s"%filename,"r")
singpos_file = open("pos-%s.vasp"%num,'w')

line = pos_file.readline()
line1 = line.strip().split()
tag = line1[0]

iswrite = False
ncount = 0
while line:
    line1 = line.strip().split()
    if len(line1) >= 2:
        if line1[1] == num:
            iswrite = True
            line1 = ['xxxxxooooo']
            ncount +=1

    if ncount == 1 and line1[0] == tag:
        iswrite = False

    if iswrite:
        singpos_file.write(line)
        
    line = pos_file.readline()

pos_file.close()
singpos_file.close()

print

print "Picked str. No.: %s from %s\n"%(num,filename)

os.system("aconvasp --std_conv <pos-%s.vasp > stdconv-%s.vasp"%(num,num))
os.system("aconvasp --std_prim <stdconv-%s.vasp > stdprim-%s.vasp\n"%(num,num))

print "Standard conventional cell is written in stdconv-%s.vasp"%num
print "Standard primitive cell is written in stdprim-%s.vasp\n"%num


POS = read_vasp("stdconv-%s.vasp"%num)

sg=spglib.get_spacegroup(POS,symprec=acc)
m=re.match('([A-Z].*\\b)\s*\(([0-9]*)\)',sg)
print "Reanalyzing symmetry using accuracy: %s A"%str(acc)
print sg
print

print "The special K points and kpath are \nwritten in k.out file according to:\n"

sgn = int(m.group(2))

a,b,c,al,be,ga = list(cell.cell_to_cellpar(POS.get_cell()))

#print sgn,a,b,c,al,be,ga

ral,rbe,rga = math.radians(al),math.radians(be),math.radians(ga)
#print ral,rbe,rga

kal = math.degrees(math.acos((math.cos(rbe)*math.cos(rga)-math.cos(ral))/(math.sin(rbe)*math.sin(rga))))
kbe = math.degrees(math.acos((math.cos(ral)*math.cos(rga)-math.cos(rbe))/(math.sin(ral)*math.sin(rga))))
kga = math.degrees(math.acos((math.cos(ral)*math.cos(rbe)-math.cos(rga))/(math.sin(ral)*math.sin(rbe))))

rkal,rkbe,rkga = math.radians(kal),math.radians(kbe),math.radians(kga)
ninterval = 21

#Triclinic:
aP = [1,2]
#monoclinic:
mS = [5,8,9,12,15]
mP = [3,4,6,7,10,11,13,14]
#Orthorhombic:
oP = [16,17,18,19,25,26,27,28,29,30,31,32,33,34,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62]
oF = [22,42,43,69,70]
oI = [23,24,44,45,46,71,72,73,74]
oS = [20,21,35,36,37,63,64,65,66,67,68,38,39,40,41]
#Hexagonal:
hP = [168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,143,144,145,147,149,150,151,152,153,154,156,157,158,159,162,163,164,165]
#Rhombohedral:
hR = [146,148,155,160,161,166,167]
#Tetragonal:
tP = [75,76,77,78,81,83,84,85,86,89,90,91,92,93,94,95,96,99,100,101,102,103,104,105,106,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138]
tI = [79,80,82,87,88,97,98,107,108,109,110,119,120,121,122,139,140,141,142]
#Cubic
cP = [195,198,200,201,205,207,208,212,213,215,218,221,222,223,224]
cF = [196,202,203,209,210,216,219,225,226,227,228]
cI = [197,199,204,206,211,214,217,220,229,230]

koutfile = open("k.out","w")

G = "0.0 0.0 0.0"
def genkp(skpts,pwscfskpts,skptssymbol):

	koutfile.write("\n")
	koutfile.write("VASP and PHON format special k-points (wrt reciprocal lattice):")
	koutfile.write("\n")
	koutfile.write("QI =")
	for i in range(len(skpts)-1):
		koutfile.write("  " + skpts[i])
		
	koutfile.write("\n")
	koutfile.write("QF =")
	
	i = 1
	while i <= len(skpts)-1:
		koutfile.write("  " + skpts[i])
		i += 1
					
	koutfile.write("\n")
	def equalsp(a,b,outfile=koutfile):
		s1 = np.linspace(a[0],b[0],ninterval)
		s2 = np.linspace(a[1],b[1],ninterval)
		s3 = np.linspace(a[2],b[2],ninterval)
		for i in range(len(s1)-1):
			outfile.write("%9.5f %9.5f %9.5f" %(s1[i],s2[i],s3[i]) +"\n")
	
	def dist(k1,k2):
		d = ((k1[0]-k2[0])**2+(k1[1]-k2[1])**2+(k1[2]-k2[2])**2)**0.5
		return d

	koutfile.write("\n")
	koutfile.write("="*10+"PWSCF input format"+"="*10+"\n")
	koutfile.write("&SYSTEM\n")
	koutfile.write("ibrav=0,\ncelldm(1)=1.8897268777743552\n")
	koutfile.write("/\n")
	koutfile.write("CELL_PARAMETERS\n")
	koutfile.write("%9.5f %9.5f %9.5f" %(array(abc)[0][0],array(abc)[0][1],array(abc)[0][2]) + "\n")
	koutfile.write("%9.5f %9.5f %9.5f" %(array(abc)[1][0],array(abc)[1][1],array(abc)[1][2]) + "\n")
	koutfile.write("%9.5f %9.5f %9.5f" %(array(abc)[2][0],array(abc)[2][1],array(abc)[2][2]) + "\n")
	koutfile.write("="*10+"PWSCF input format"+"="*10+"\n")
#	koutfile.write("\n")
#	koutfile.write("Reciprocal vectors:" + "\n")
#	rabc = (abc.I).T
#	koutfile.write("%9.5f %9.5f %9.5f" %(array(rabc)[0][0],array(rabc)[0][1],array(rabc)[0][2]) + "\n")
#	koutfile.write("%9.5f %9.5f %9.5f" %(array(rabc)[1][0],array(rabc)[1][1],array(rabc)[1][2]) + "\n")
#	koutfile.write("%9.5f %9.5f %9.5f" %(array(rabc)[2][0],array(rabc)[2][1],array(rabc)[2][2]) + "\n")
	koutfile.write("\n")
	koutfile.write("PWSCF format special k-points (Cartesian in reciprocal space):")
	koutfile.write("\n")
	dkpts = 0
	for i in range(len(pwscfskpts)-1):
		koutfile.write("%9.5f %9.5f %9.5f   @%9.5f  %s" %(pwscfskpts[i][0][0],pwscfskpts[i][0][1],pwscfskpts[i][0][2],dkpts,skptssymbol[i]) + "\n")
		dkpts += dist(pwscfskpts[i][0],pwscfskpts[i+1][0])

	koutfile.write("%9.5f %9.5f %9.5f   @%9.5f  %s" %(pwscfskpts[len(pwscfskpts)-1][0][0],pwscfskpts[len(pwscfskpts)-1][0][1],pwscfskpts[len(pwscfskpts)-1][0][2],dkpts,skptssymbol[len(pwscfskpts)-1]) + "\n")
	koutfile.write("\n")
	koutfile.write("PWSCF format k-path (Cartesian in reciprocal space):")
	koutfile.write("\n")
	koutfile.write("%.0f" %((len(skpts)-1)*(ninterval-1)+1) +"\n")
	for i in range(len(pwscfskpts)-1):
		equalsp(pwscfskpts[i][0],pwscfskpts[i+1][0])

	koutfile.write("%9.5f %9.5f %9.5f" %(float(pwscfskpts[len(skpts)-1][0][0]), float(pwscfskpts[len(skpts)-1][0][1]),float(pwscfskpts[len(skpts)-1][0][2])))
	
	
if sgn in aP:
	abc = matrix([[a,0,0],[b*math.cos(rga),b*math.sin(rga),0],[c*math.cos(rbe),c*(math.cos(ral)-math.cos(rbe)*math.cos(rga))/math.sin(rga),c/math.sin(rga)*(math.sin(rga)**2-math.cos(ral)**0.5-math.cos(rbe)**2+2*math.cos(ral)*math.cos(rbe)*math.cos(rga))**0.5]])
	if kal > 90 and kbe > 90 and kga > 90 and	kga == min(kal,kbe,kga):
		X = "0.5 0.0 0.0"
		Y = "0.0 0.5 0.0"
		koutfile.write("X-Gamma-Y" + "\n")
		skpts = [X,G,Y]
		skptssymbol = ["X","G","Y"]
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "aP1a"

	elif kal > 90 and kbe > 90 and kga == 90:
		X = "0.5 0.0 0.0"
		Y = "0.0 0.5 0.0"
		koutfile.write("X-Gamma-Y" + "\n")
		skpts = [X,G,Y]
		skptssymbol = ["X","G","Y"]
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "aP2a"

	elif kal < 90 and kbe < 90 and kga < 90 and kga == max(kal,kbe,kga):
		X = "0.0 -0.5 0.0"
		Y = "0.5 0.0 0.0"
		koutfile.write("X-Gamma-Y" + "\n")
		skpts = [X,G,Y]
		skptssymbol = ["X","G","Y"]
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "aP1b"

	elif kal < 90 and kbe < 90 and kga == 90:
		X = "0.0 -0.5 0.0"
		Y = "0.5 0.0 0.0"
		koutfile.write("X-Gamma-Y" + "\n")
		skpts = [X,G,Y]
		skptssymbol = ["X","G","Y"]
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "aP2b"
		
	else:
		print "aP1a"

elif sgn in mP or sgn in mS:
	ll = [a,b,c]
	ll.sort()
	aa = ll[0]
	bb = ll[1]
	cc = ll[2]
	degs = [al,be,ga]
	degs.sort()
	alpha = degs[0]
	beta = degs[1]
	gamma = degs[2]

	if gamma >= 90:
		alpha = 180-gamma
		gamma = 90
		oo = bb
		bb = aa
		aa = oo
		
	ralpha = math.radians(alpha)
	rbeta = math.radians(beta)
	rgamma = math.radians(gamma)
	
	if sgn in mP:
		eta = (1-bb*math.cos(ralpha)/cc)/(2*(math.sin(ralpha))**2)
		nu = 0.5-eta*cc*math.cos(ralpha)/bb
		koutfile.write("NOTE: a,b,c must be in ascending order,alpha <90,beta=gamma=90!!" + "\n")
		koutfile.write("a= %8.6f, b= %8.6f, c= %8.6f, alpha= %6.4f,	beta= %6.4f, gamma= %6.4f" %(aa,bb,cc,alpha,beta,gamma) + "\n\n")
		Y = "0.0 0.0 0.5"
		H = "0.0 %9.5f %9.5f"%(eta,1.0-nu)
		C = "0.0 0.5 0.5"
		E = "0.5 0.5 0.5"
		M1= "0.5 %9.5f %9.5f"%(1.0-eta,nu)
		A = "0.5 0.5 0.0 "
		X = "0.0 0.5 0.0"
		H1= "0.0 %9.5f %9.5f"%(1.0-eta,nu)
		koutfile.write("Gamma-Y-H-C-E-M1-A-X-H1" + "\n")
		skpts = [G,Y,H,C,E,M1,A,X,H1]
		skptssymbol = ["G","Y","H","C","E","M1","A","X","H1"]
		abc = matrix([[aa,0,0],[0,bb,0],[0,cc*math.cos(ralpha),cc*math.sin(ralpha)]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "mP"

	elif sgn in mS:
		koutfile.write("NOTE: a,b,c must be in ascending order,alpha <90,beta=gamma=90!!" + "\n")
		koutfile.write("a= %8.6f, b= %8.6f, c= %8.6f, alpha= %9.5f, beta= %9.5f, gamma= %9.5f" %(aa,bb,cc,alpha,beta,gamma) + "\n\n")
		abc = matrix([[aa/2,bb/2,0],[-aa/2,bb/2,0],[0,cc*math.cos(ralpha),cc*math.sin(ralpha)]])

		if kga > 90:
			zeta = (2.0-bb*math.cos(ralpha)/cc)/(4*(math.sin(ralpha))**2)
			eta = 0.5+2*zeta*cc*math.cos(ralpha)/bb
			psi = 0.75-aa**2/(4*bb**2*(math.sin(ralpha))**2)
			phi = psi+(0.75-psi)*bb**2*math.cos(ralpha)/cc
			Y = "0.5 0.5 0.0"
			F = "%9.5f %9.5f %9.5f"%(1.0-zeta,1.0-zeta,1.0-zeta)
			L = "0.5 0.5 0.5"
			I = "%9.5f %9.5f 0.5"%(phi,1.0-phi)
			koutfile.write("Gamma-Y-F-L-I" + "\n")
			skpts = [G,Y,F,L,I]
			skptssymbol = ["G","Y","F","L","I"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print  "mS1"

		elif kga == 90:
			zeta = (2.0-bb*math.cos(ralpha)/cc)/(4*(math.sin(ralpha))**2)
			eta = 0.5+2*zeta*cc*math.cos(ralpha)/bb
			psi = 0.75-aa**2/(4*bb**2*(math.sin(ralpha))**2)
			phi = psi+(0.75-psi)*bb**2*math.cos(ralpha)/cc
			
			Y = "0.5 0.5 0.0"
			F = "%9.5f %9.5f %9.5f"%(1.0-zeta,1.0-zeta,1.0-zeta)
			L = "0.5 0.5 0.5"
			I = "%9.5f %9.5f 0.5"%(phi,1.0-phi)
			koutfile.write("Gamma-Y-F-L-I" + "\n")
			skpts = [G,Y,F,L,I]
			skptssymbol = ["G","Y","F","L","I"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print  "mS2"

		elif kga < 90 and bb*math.cos(ralpha)/cc+bb**2*(math.sin(ralpha))**2/aa**2 < 1:
			mu = (1+bb*2/aa*2)/4.0
			delta = bb*cc*math.cos(ralpha)/(2*aa**2)
			zeta = mu-0.25+(1-bb*math.cos(ralpha)/cc)/(4*(math.sin(ralpha))**2)
			eta = 0.5+2*zeta*cc*math.cos(ralpha)/bb
			phi = 1+zeta-2*mu
			psi = eta-2*delta
			Y = "%9.5f %9.5f %9.5f"%(mu,mu,delta)
			F = "%9.5f %9.5f %9.5f"%(1.0-phi,1.0-phi,1.0-psi)
			H = " %9.5f %9.5f %9.5f"%(zeta,zeta,eta)
			Z = "0.0 0.0 0.5"
			I = "0.5 -0.5 0.5"
			F1 = "%9.5f %9.5f %9.5f"%(phi,phi-1.0,psi)
			koutfile.write("Gamma-Y-F-H-Z-I-F1" + "\n")
			skpts = [G,Y,F,H,Z,I,F1]
			skptssymbol = ["G","Y","F","H","Z","I","F1"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "mS3"

		elif kga < 90 and bb*math.cos(ralpha)/cc+bb**2*(math.sin(ralpha))**2/aa**2 == 1:
			mu = (1+bb*2/aa*2)/4.0
			delta = bb*cc*math.cos(ralpha)/(2*aa**2)
			zeta = mu-0.25+(1-bb*math.cos(ralpha)/cc)/(4*(math.sin(ralpha))**2)
			eta = 0.5+2*zeta*cc*math.cos(ralpha)/bb
			phi = 1+zeta-2*mu
			psi = eta-2*delta
			Y = "%9.5f %9.5f %9.5f"%(mu,mu,delta)
			F = "%9.5f %9.5f %9.5f"%(1.0-phi,1.0-phi,1.0-psi)
			H = "%9.5f %9.5f %9.5f"%(zeta,zeta,eta)
			Z = "0.0 0.0 0.5"
			I = "0.5 -0.5 0.5"
			koutfile.write("Gamma-Y-F-H-Z-I" + "\n")
			skpts = [G,Y,F,H,Z,I]
			skptssymbol = ["G","Y","F","H","Z","I"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "mS4"

		elif kga < 90 and bb*math.cos(ralpha)/cc+bb**2*(math.sin(ralpha))**2/aa**2 > 1:
			zeta = (bb**2/aa**2+(1-bb*math.cos(ralpha)/cc)/(math.sin(ralpha))**2)/4.0
			eta = 0.5+2.0*zeta*cc*math.cos(ralpha)/bb
			mu = eta/2.0+bb**2/4.0/aa*2-bb*cc*math.cos(ralpha)/2.0/aa*2
			nu = 2*mu-zeta
			omega = (4.0*nu-1.0-bb*2*(math.sin(ralpha))**2/aa**2)*cc/(2*bb*math.cos(ralpha))
			delta = zeta*cc*math.cos(ralpha)/bb+omega/2.0-0.25
			rho = 1.0-zeta*aa**2/bb**2
			Y = "%9.5f %9.5f %9.5f"%(mu,mu,delta)
			F = "%9.5f %9.5f %9.5f"%(nu,nu,omega)
			L = "0.5 0.5 0.5"
			I = "%9.5f %9.5f 0.5"%(rho,1.0-rho)
			koutfile.write("Gamma-Y-F-L-I" + "\n")
			skpts = [G,Y,F,L,I]
			skptssymbol = ["G","Y","F","L","I"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "mS5"
		

elif sgn in oP or sgn in oF or sgn in oI or sgn in oS:
	ll = [a,b,c]
	ll.sort()
	aa = ll[0]
	bb = ll[1]
	cc = ll[2]

	if sgn in oS:
		zeta = (1+a**2/b**2)/4.0
		X = "%9.5f  %9.5f 0.0"%(zeta,zeta)
		S = "0.0 0.5 0.0"
		R = "0.0 0.5 0.5"
		A = "%9.5f  %9.5f  0.5"%(zeta,zeta)
		Z = "0.0 0.0 0.5"
		Y = "-0.5 0.5 0.0"
		X1= "%9.5f  %9.5f 0.0"%(-zeta,1.0-zeta)
		A1= "%9.5f  %9.5f 0.5"%(-zeta,1.0-zeta)
		T = "-0.5 0.5 0.5"
		Y = "-0.5 0.5 0.0"
		koutfile.write("Gamma-X-S-R-A-Z-Gamma-Y-X1-A1-T-Y" + "\n")
		skpts = [G,X,S,R,A,Z,G,Y,X1,A1,T,Y]
		skptssymbol = ["G","X","S","R","A","Z","G","Y","X1","A1","T","Y"]
		abc = matrix([[aa/2,-bb/2,0],[aa/2,bb/2,0],[0,0,cc]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "oS"


	if sgn in oP:
		koutfile.write("NOTE: a,b,c must be in ascending order!!" + "\n")
		koutfile.write("a= %8.6f, b= %8.6f, c= %8.6f, alpha= %6.4f, beta= %6.4f, gamma= %6.4f" %(aa,bb,cc,al,be,ga) + "\n\n")
		abc = matrix([[aa,0,0],[0,bb,0],[0,0,cc]])
		X = "0.5 0.0 0.0"
		S = "0.5 0.5 0.0"
		Y = "0.0 0.5 0.0"
		Z = "0.0 0.0 0.5"
		U = "0.5 0.0 0.5"
		R = "0.5 0.5 0.5"
		T = "0.0 0.5 0.5"
		koutfile.write("Gamma-X-S-Y-Gamma-Z-U-R-T-Z" + "\n")
		skpts = [G,X,S,Y,G,Z,U,R,T,Z]
		skptssymbol = ["G","X","S","Y","G","Z","U","R","T","Z"]
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print  "oP"

	elif sgn in oF:
		koutfile.write("NOTE: a,b,c must be in ascending order!!" + "\n")
		koutfile.write("a= %8.6f, b= %8.6f, c= %8.6f, alpha= %6.4f,	beta= %6.4f, gamma= %6.4f" %(aa,bb,cc,al,be,ga) + "\n\n")
		abc = matrix([[0,bb/2,cc/2],[aa/2,0,cc/2],[aa/2,bb/2,0]])
		if 1/aa**2 > 1/bb**2 + 1/cc**2:
			zeta =(1+aa**2/bb**2-aa**2/cc**2)/4.0
			eta = (1+aa**2/bb**2+aa**2/cc**2)/4.0
			Y = "0.5 0.0 0.5"
			T = "1.0 0.5 0.5"
			Z = "0.5 0.5 0.0"
			X = "0.0 %9.5f %9.5f"%(eta,eta)
			A1= "0.5 %9.5f %9.5f"%(0.5-zeta,1.0-zeta)
			koutfile.write("Gamma-Y-T-Z-Gamma-X-A1-Y" + "\n")
			skpts = [G,Y,T,Z,G,X,A1,Y]
			skptssymbol = ["G","Y","T","Z","G","X","A1","Y"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "oF1"

		elif 1/aa**2 < 1/bb**2 + 1/cc**2:
			eta = (1+aa**2/bb**2-aa**2/cc**2)/4.0
			phi = (1+cc**2/bb**2-cc**2/aa**2)/4.0
			delta = (1+bb**2/aa**2-bb**2/cc**2)/4.0
			Y = "0.5 0.0 0.5"
			C = "0.5 %9.5f %9.5f"%(0.5-eta,1.0-eta)
			D = "%9.5f 0.5 %9.5f"%(0.5-delta,1.0-delta)
			X = "0.0 0.5 0.5"
			Z = "0.5 0.5 0.0"
			D1 = "%9.5f 0.5 %9.5f"%(0.5+delta,delta)
			H = "%9.5f %9.5f 0.5"%(1.0-phi,0.5-phi)
			C = "0.5 %9.5f %9.5f"%(0.5-eta,1.0-eta)
			koutfile.write("Gamma-Y-C-D-X-Gamma-Z-D1-H-C" + "\n")
			skpts = [G,Y,C,D,X,G,Z,D1,H,C]
			skptssymbol = ["G","Y","C","D","X","G","Z","D1","H","C"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "oF2"

		elif 1/aa**2 == 1/bb**2 + 1/cc**2:
			zeta == (1+aa**2/bb**2-aa**2/cc**2)/4.0
			eta == (1+aa**2/bb**2+aa**2/cc**2)/4.0
			Y = "0.5 0.0 0.5"
			T = "1.0 0.5 0.5"
			Z = "0.5 0.5 0.0"
			X = "0.0 %9.5f %9.5f"%(eta,eta)
			A1= "0.5 %9.5f %9.5f"%(0.5-zeta,1.0-zeta)
			koutfile.write("Gamma-Y-T-Z-Gamma-X-A1-Y" + "\n")
			skpts = [G,Y,T,Z,G,X,A1,Y]
			skptssymbol = ["G","Y","T","Z","G","X","A1","Y"]
			spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
			pwscfskpts = [i for i in array(spp)]
			genkp(skpts,pwscfskpts,skptssymbol)
			print "oF3"

	elif sgn in oI:
		koutfile.write("NOTE: a,b,c must be in ascending order!!" + "\n")
		koutfile.write("a= %8.6f, b= %8.6f, c= %8.6f, alpha= %6.4f,	beta= %6.4f, gamma= %6.4f" %(aa,bb,cc,al,be,ga) + "\n\n")
		zeta = (1+aa**2/cc**2)/4.0
		eta = (1+bb**2/cc**2)/4.0
		delta = (bb**2-aa**2)/(4.0*cc**2)
		mu = (aa**2+bb**2)/(4.0*cc**2)
		X = "%9.5f  %9.5f  %9.5f"%(-zeta,zeta,zeta)
		L = "%9.5f  %9.5f  %9.5f"%(-mu,mu,0.5-delta)
		T = "0.0 0.0 0.5"
		W = "0.25 0.25 0.25"
		R = "0.0 0.5 0.0"
		X1= "%9.5f  %9.5f  %9.5f"%(zeta,1.0-zeta,-zeta)
		Z = "0.5 0.5 -0.5"
		Y = "%9.5f  %9.5f  %9.5f"%(eta,-eta,eta)
		S = "0.5 0.0 0.0"
		koutfile.write("Gamma-X-L-T-W-R-X1-Z-Gamma-Y-S-W" + "\n")
		skpts = [G,X,L,T,W,R,X1,Z,G,Y,S,W]
		skptssymbol = ["G","X","L","T","W","R","X1","Z","G","Y","S","W"]
		abc = matrix([[-aa/2,bb/2,cc/2],[aa/2,-bb/2,cc/2],[aa/2,bb/2,-cc/2]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "oI"


elif sgn in hP:
	M = "0.5 0.0 0.0"
	K = "0.33333 0.33333 0.0"
	A = "0.0 0.0 0.5"
	L = "0.5 0.0 0.5"
	H = "0.33333 0.33333 0.5"
	koutfile.write("Gamma-M-K-Gamma-A-L-H-A" + "\n")
	skpts = [G,M,K,G,A,L,H,A]
	skptssymbol = ["G","M","K","G","A","L","H","A"]
	abc = matrix([[a/2,-a*3**0.5/2,0],[a/2,a*3**0.5/2,0],[0,0,c]])
	spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
	pwscfskpts = [i for i in array(spp)]
	genkp(skpts,pwscfskpts,skptssymbol)
	print "hP"

elif sgn in hR:
	if al < 90:
		eta = (1+4*math.cos(ral))/(2.0+4*math.cos(ral))
		nu = 3.0/4.0-eta/2
		L = "0.5 0.0 0.0"
		B1= "0.5 %9.5f %9.5f"%(1.0-eta,eta-1.0)
		koutfile.write("Gamma-L-B1" + "\n")
		skpts = [G,L,B1]
		skptssymbol = ["G","L","B1"]
		abc = matrix([[a*math.cos(ral/2),-a*math.sin(ral/2),0],[a*math.cos(ral/2),a*math.sin(ral/2),0],[a*math.cos(ral)/math.cos(ral/2),0,a*(1-(math.cos(ral)/math.cos(ral/2))**2)**0.5]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "hR1"

	elif al > 90:
		eta = 1/(2*(math.tan(ral/2.0))**2)
		nu = 3.0/4.0-eta/2.0
		P = "%9.5f %9.5f %9.5f"%(1.0-nu,-nu,1.0-nu)
		Z = "0.5 -0.5 0.5"
		Q = "%9.5f %9.5f %9.5f"%(eta,eta,eta)
		F = "0.5 -0.5 0.0"
		P1="%9.5f %9.5f %9.5f"%(nu,nu-1.0,nu-1.0)
		Q1="%9.5f %9.5f %9.5f"%(1.0-eta,-eta,-eta)
		L="0.5 0.0 0.0"
		koutfile.write("Gamma-P-Z-Q-Gamma-F-P1-Q1-L-Z" + "\n")
		skpts = [G,P,Z,Q,G,F,P1,Q1,L,Z]
		skptssymbol = ["G","P","Z","Q","G","F","P1","Q1","L","Z"]
		abc = matrix([[a*math.cos(ral/2),-a*math.sin(ral/2),0],[a*math.cos(ral/2),a*math.sin(ral/2),0],[a*math.cos(ral)/math.cos(ral/2),0,a*(1-(math.cos(ral)/math.cos(ral/2))**2)**0.5]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print "hR2"

elif sgn in tP:
	X = "0.0 0.5 0.0"
	M = "0.5 0.5 0.0"
	Z = "0.0 0.0 0.5"
	R = "0.0 0.5 0.5"
	A = "0.5 0.5 0.5"
	koutfile.write("Gamma-X-M-G-Z-R-A-Z" + "\n")
	skpts = [G,X,M,G,Z,R,A,Z]
	skptssymbol = ["G","X","M","G","Z","R","A","Z"]
	abc = matrix([[a,0,0],[0,a,0],[0,0,c]])
	spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
	pwscfskpts = [i for i in array(spp)]
	genkp(skpts,pwscfskpts,skptssymbol)
	print "tP"

elif sgn in tI:
	if c < a:
		eta = (1.0 + c**2/a**2)/4.0
		X = "0.0 0.0 0.5"
		M = "-0.5 0.5 0.5"
		Z = "%9.5f %9.5f %9.5f"%(eta,eta,-eta)
		P = "0.25 0.25 0.25"
		N = "0.0 0.5 0.0"
		Z1= "%9.5f %9.5f %9.5f"%(-eta,1.0-eta,eta)
		koutfile.write("Gamma-X-M-Gamma-Z-P-N-Z1-M" + "\n")
		skpts = [G,X,M,G,Z,P,N,Z1,M]
		skptssymbol = ["G","X","M","G","Z","P","N","Z1","M"]
		abc = matrix([[-a/2,a/2,c/2],[a/2,-a/2,c/2],[a/2,a/2,-c/2]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print  "tI1"

	elif c > a:
		eta = (1.0 + a**2/c**2)/4.0
		zeta = a**2/(2*c**2)
		X = "0.0 0.0 0.5"
		Y = "%9.5f %9.5f 0.5"%(-zeta,zeta)
		S = "%9.5f %9.5f %9.5f"%(-eta,eta,eta)
		Z = "0.5 0.5 -0.5"
		S1= "%9.5f %9.5f %9.5f"%(eta,1-eta,-eta)
		N = "0.0 0.5 0.0"
		P = "0.25 0.25 0.25"
		Y1= "0.5 0.5 %9.5f"%(-zeta)
		koutfile.write("Gamma-X-Y-Sigma-Gamma-Z-Sigma1-N-P-Y1-Z" + "\n")
		skpts = [G,X,Y,S,G,Z,S1,N,P,Y1,Z]
		skptssymbol = ["G","X","Y","S","G","Z","S1","N","P","Y1","Z"]
		abc = matrix([[-a/2,a/2,c/2],[a/2,-a/2,c/2],[a/2,a/2,-c/2]])
		spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
		pwscfskpts = [i for i in array(spp)]
		genkp(skpts,pwscfskpts,skptssymbol)
		print  "tI2"

elif sgn in cP:
	X = "0.0 0.5 0.0"
	M = "0.5 0.5 0.0"
	R = "0.5 0.5 0.5"
	koutfile.write("G-X-M-G-R-X" + "\n")
	skpts = [G,X,M,G,R,X]
	skptssymbol = ["G","X","M","G","R","X"]
	abc = matrix([[a,0,0],[0,a,0],[0,0,a]])
	spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
	pwscfskpts = [i for i in array(spp)]
	genkp(skpts,pwscfskpts,skptssymbol)
	print "cP"

elif sgn in cF:
	X = "0.500 0.000 0.500"
	W = "0.500 0.250 0.750"
	K = "0.375 0.375 0.750"
	L = "0.500 0.500 0.500"
	U = "0.625 0.250 0.625"
	koutfile.write("G-X-W-K-G-L-U-W-L-K" + "\n")
	skpts = [G,X,W,K,G,L,U,W,L,K]
	skptssymbol = ["G","X","W","K","G","L","U","W","L","K"]
	abc = matrix([[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0]])
	spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
	pwscfskpts = [i for i in array(spp)]
	genkp(skpts,pwscfskpts,skptssymbol)
	print "cF"

elif sgn in cI:
	H = "0.50 -0.50 0.50"
	N = "0.00  0.00 0.50"
	P = "0.25  0.25 0.25"
	koutfile.write("Gamma-H-N-Gamma-P-H" + "\n")
	skpts = [G,H,P,G,N]
	skptssymbol = ["G","H","P","G","N"]
	abc = matrix([[-a/2,a/2,a/2],[a/2,-a/2,a/2],[a/2,a/2,-a/2]])
	spp = [i*((abc.I).T) for i in [matrix(j) for j in skpts]]
	pwscfskpts = [i for i in array(spp)]
	genkp(skpts,pwscfskpts,skptssymbol)
	print "cI"

else:
	print "BAD NEWS!! Space Group Number is wrong."

print "\nsymmetry."
