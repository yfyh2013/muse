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

from math import sin,cos,pi,sqrt
from muse.Symmetry.Const import *
from muse.Symmetry.SG import *
import random
from muse.Genstrs import GENCELL,Gen_sg_coord,Generate_Supercell
from muse.Readwrite.ReadInput import indict

class GenPC(GeometryObject):
    def __init__(self,sg_num,cellp,occ,ineqsites):
        self.a = cellp[0]
        self.b = cellp[1]
        self.c = cellp[2]
        self.alpha = cellp[3]
        self.beta = cellp[4]
        self.gamma = cellp[5]
        self.spacegroupnr = sg_num
        self.HMSymbol = ""
        self.ineqsites = ineqsites
        self.atomdata = []
        self.atomset = set([])
        self.occupations = occ
        #self.latticevectors = None
        self.lengthscale = 1
        self.symops = set([])
        # Lattice vector choices
        self.primcell = False
        self.rhomb2hex = False
        self.rhombohedral = False
        self.spacegroupsetting = ""
        self.coordepsilon = 0.0002
        
    def cry_system(self,spacegroupnr):
        # Determine crystal system
        if 0 < spacegroupnr <= 2:
            return "triclinic"
        elif 2 < spacegroupnr <= 15:
            return "monoclinic"
        elif 15 < spacegroupnr <= 74:
            return "orthorhombic"
        elif 74 < spacegroupnr <= 142:
            return "tetragonal"
        elif 142 < spacegroupnr <= 167:
            return "trigonal"
        elif 167 < spacegroupnr <= 194:
            return "hexagonal"
        elif 194 < spacegroupnr <= 230:
            return "cubic"
        else:
            return "unknown"

    def crystal_system(self):
        return self.cry_system(self.spacegroupnr)

    def conventional_latticevectors(self):
        # Set up Bravais lattice vectors of the conventional cell
        self.coa = self.c / self.a
        self.boa = self.b / self.a
        alphar = self.alpha*pi/180
        betar  = self.beta*pi/180
        gammar = self.gamma*pi/180
        if self.crystal_system() == 'cubic':
            latticevectors = LatticeMatrix([[one, zero, zero],
                                            [zero, one, zero],
                                            [zero, zero, one]])
        elif self.crystal_system() == 'hexagonal':
            latticevectors = LatticeMatrix([[sin(gammar), cos(gammar), zero],
                                            [zero, one, zero],
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'tetragonal' or self.crystal_system() == 'orthorhombic':
            latticevectors = LatticeMatrix([[one, zero, zero], 
                                            [zero, self.boa, zero], 
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'monoclinic':
            latticevectors = LatticeMatrix([[one, zero, zero],
                                            [zero, self.boa, zero],
                                            [self.coa*cos(betar), zero, self.coa*sin(betar)]])
        elif self.crystal_system() == 'trigonal':
            # Hexagonal cell taken as conventional
            if not abs(self.gamma-120) < self.coordepsilon:
                gammar = 120*pi/180
            latticevectors = LatticeMatrix([[sin(gammar), cos(gammar), zero],
                                            [zero, one, zero],
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'triclinic' or self.crystal_system() == 'monoclinic' or self.crystal_system() == 'unknown':
            angfac1 = (cos(alphar) - cos(betar)*cos(gammar))/sin(gammar)
            angfac2 = sqrt(sin(gammar)**2 - cos(betar)**2 - cos(alphar)**2 
                       + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)
            latticevectors = LatticeMatrix([[one, zero, zero], 
                                            [self.boa*cos(gammar), self.boa*sin(gammar), zero], 
                                            [self.coa*cos(betar), self.coa*angfac1, self.coa*angfac2]])
        else:
            raise SymmetryError("No support for "+self.crystal_system()+" crystal systems.")
        return latticevectors
        
    def getCrystalStructure(self,IfReduce=True):
        self.primcell = IfReduce
#        print self.spacegroupnr
#        print "="*100

        if not 0 < self.spacegroupnr <= 230:
            try:
                self.spacegroupnr = Hall2Number[self.HallSymbol]
            except:
                pass
        if self.HMSymbol == "":
            try:
                self.HMSymbol == NumberHM[self.HallSymbol]
            except:
                pass

        if self.spacegroupsetting == "":
            try:
                self.spacegroupsetting = self.HallSymbol.lstrip("-")[0]
            except:
                try:
                    self.spacegroupsetting = self.HMSymbol[0]
                except:
                    pass
        if self.spacegroupsetting == "":
            try:
                self.spacegroupsetting = SGnrtoHM[str(self.spacegroupnr)][0]
            except:
                pass
        # Sanity test of lattice parameters
        if self.a!=0 and self.b!=0 and self.c!=0 and self.alpha!=0 and self.beta!=0 and self.gamma!=0:
            if not 0 < self.spacegroupnr < 231:
                if len(self.symops) >= 1:
                    if self.primcell == True:
                        if self.spacegroupsetting == 'P':
                            self.spacegroupnr = 0
                        else:
                            raise SymmetryError("Insufficient symmetry information to reduce to primitive cell"+\
                                                " (need space group number or Hermann-Mauguin symbol).\n"+\
                                                "  Run with --no-reduce to generate cell in the conventional setting.")
                    else:
                        self.spacegroupnr = 0
        else:
            raise CellError("No crystallographic parameter may be zero.")

        # If no symmetry operations are not set, get internally stored.
        if not self.symops:
            eqsites = SymmetryOp[self.spacegroupnr]
            # Define the set of space group operations.
            self.symops = set([])
            for site in eqsites:
                self.symops.add(SymmetryOperation(site))

        self.latticevectors = self.conventional_latticevectors()
        self.lengthscale = self.a

        if self.primcell:
            if self.spacegroupsetting == 'I':
                # Body centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,half])]
                if self.crystal_system() == 'cubic':
                    self.lattrans = LatticeMatrix([[-half, half, half],
                                                   [half, -half, half],
                                                   [half, half, -half]])
                else:
                    self.lattrans = LatticeMatrix([[one, zero, zero],
                                                   [zero, one, zero],
                                                   [half, half, half]])
            elif self.spacegroupsetting == 'F':
                # face centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero]),
                                  LatticeVector([half,zero,half]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = LatticeMatrix([[half, half, zero],
                                               [half, zero, half],
                                               [zero, half, half]])
            elif self.spacegroupsetting == 'A':
                # A-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = LatticeMatrix([[one, zero, zero],
                                               [zero, half, -half],
                                               [zero, half, half]])
            elif self.spacegroupsetting == 'B':
                # B-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,zero,half])]
                self.lattrans = LatticeMatrix([[half, zero, -half],
                                               [zero, one, zero],
                                               [half, zero, half]])
            elif self.spacegroupsetting == 'C':
                # C-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero])]
                self.lattrans = LatticeMatrix([[half, -half, zero],
                                               [half, half, zero],
                                               [zero, zero, one]])
            else:
                self.transvecs = [LatticeVector([zero,zero,zero])]
                self.lattrans = LatticeMatrix([[1, 0, 0],
                                               [0, 1, 0],
                                               [0, 0, 1]])
            # Transform to primitive cell
            tmp = []
            for i in range(3):
                tmp.append(mvmult3(self.latticevectors,self.lattrans[i]))
            self.latticevectors = LatticeMatrix(tmp)
            # Improve precision again...
            
            for i in range(3):
                for j in range(3):
                    self.latticevectors[i][j] = improveprecision(self.latticevectors[i][j],self.coordepsilon)
        else:
            # If no reduction is to be done
            self.transvecs = [LatticeVector([zero,zero,zero])]
            self.lattrans = LatticeMatrix([[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])
                                          
        # Find inverse lattice transformation matrix
        invlattrans = LatticeMatrix(minv3(self.lattrans))
        

        if self.primcell:
            if len(self.transvecs) > 1:
                redundant = set([])
                for op1 in self.symops:
                    for op2 in self.symops:
                        for vec in self.transvecs:
                            if op1.translation+vec == op2.translation:
                                if op1.rotation == op2.rotation:
                                    if op1.translation.length() < op2.translation.length():
                                        redundant.add(op2)
                self.symops -= redundant

        # Space group operations to cartesian representation
        lv = self.conventional_latticevectors()
        
        for op in self.symops:
            op.rotation = lv.transform(op.rotation)
            op.rotation = op.rotation.transform(minv3(lv))
            # transform translations
            op.translation = op.translation.transform(minv3(self.lattrans))

        for i in range(len(self.ineqsites)):
            # Set up atomdata
            self.atomdata.append([])
            self.atomdata[i].append(AtomSite(position=self.ineqsites[i]))
            # Add species and occupations to atomdata
            for k,v in self.occupations[i].iteritems():
                self.atomdata[i][0].species[k] = v
                # Add charge state

            ## self.atomdata[i][0].charge = self.charges[i]
            # Determine if we have an alloy
            for element in self.occupations[i]:
                v = self.occupations[i][element]
                if abs(1-v) > occepsilon:
                    self.alloy = True

        # Work out all sites in the cell for atomdata/atomset
        for a in self.atomdata:
            for op in self.symops:
                # position expression string
                posexpr = [s for s in op.eqsite]
                for k in range(3):
                    # position expression string, replacing x,y,z with numbers
                    posexpr[k] = posexpr[k].replace('x',str(a[0].position[0]))
                    posexpr[k] = posexpr[k].replace('y',str(a[0].position[1]))
                    posexpr[k] = posexpr[k].replace('z',str(a[0].position[2]))
                position = LatticeVector([safe_matheval(pos) for pos in posexpr])
                b = AtomSite(position=position,species=a[0].species,charges=a[0].charges)
                self.atomset.add(b)
                append = True
                for site in a:
                    for vec in self.transvecs:
                        t = vec + b.position
                        if site.position == t:
                            append=False
                            break
                    if not append:
                        break
                if append:
                    a.append(b)
        # Transform positions. Note that atomdata and atomset alias the same data,
        # so we only transform once.
        for a in self.atomdata:
            for b in a:
                b.position = LatticeVector(mvmult3(invlattrans,b.position))

        spec = []
        self.initialized = True
        scl_pos = []
        for i in range(len(self.atomdata)):
            for j in range(len(self.atomdata[i])):
                spec.append(self.atomdata[i][j].spcstring())
                scl_pos.append(self.atomdata[i][j].position)
        
        cellread = []
        for ii in self.latticevectors:
            cellread.append(ii)

        primcell = Generate_Supercell.Generate_Supercell(spec, cell=cellread, scaled_positions=scl_pos)
        primcell.set_cell(primcell.get_cell()*self.a,scale_atoms=True)
        return primcell


if __name__ == '__main__':
    nu = 5
    Num_atomsInChemForm = [int(i) for i in indict['Num_atomsInChemForm']]
    oldlist = [nu*i for i in Num_atomsInChemForm]
    iffind = False
    while not iffind:
        colist = []
        while colist == []: 
            Hall = random.choice(SymmetryOp.keys())
            sg_num = self.spacegroupnr
            colist,Nspec,issuccess = Gen_sg_coord.Gen_sg_coord(sg_num,Num_atomsInChemForm)
        
        print sg_num
        ntry = 0
        while True:
            gencell = GENCELL.GENCELL(sg_num)
            cellp = gencell.Gen_sg_latt()
            
            while min(cellp[0],cellp[1],cellp[2]) < 1.0:
                for i in range(3): cellp[i] *= 10
            
            #print cellp
            
            occ = [{indict['NameOfAtoms'][i]:Nspec[i]} for i in range(len(indict['NameOfAtoms']))]
            #print occ
            ineqsites = [list(colist[i]) for i in range(len(colist))]
            #print ineqsites
            genpc = GenPC(sg_num,Hall,cellp,occ,ineqsites)
            primcell = genpc.getCrystalStructure(IfReduce=True)
            #print "Latt. Vector:"
            #print primcell.get_cell()
            #print "Pos.:"
            #print primcell.get_scaled_positions()
            #print primcell.get_chemical_symbols()
            b = set(primcell.get_chemical_symbols())
            newlist = []
            for a in b:
                newlist.append(primcell.get_chemical_symbols().count(a))
            
            ifone = []
            ifdive = []
            for i,num in enumerate(oldlist):
                ifone.append(float(oldlist[i]/newlist[i]))
                ifdive.append(oldlist[i]%newlist[i] == 0)
        
            #print ifdive
            #print len(set(ifone))
            #print ifone[0]
            #print "oldlist=",oldlist
            print "newlist=",newlist
            ntry +=1
            if ntry > 10 or (len(set(ifone)) == 1.0 and ifone[0] >= 1.0):
                if len(set(ifone)) == 1.0 and ifone[0] >= 1.0 and len(set(ifdive)) == 1.0 and ifdive[0]:
                    iffind = True
                break
