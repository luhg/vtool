#!/usr/bin/env python

from atm import *
import re
import copy

__author__ = ""
__date__ = "2013/12/19"

__version__ = "$Revision: 0.2$"


class POSCAR:
    """ create/read VASP POSCAR """
    def __init__(self, filename = None, comment = None, lattice = None,
                 select="Selective", coordinate="Cartesian"):

        # setup comment
        if comment is None:
            self.setComment("Comment line")
        else:
            self.setComment(comment)
        self.setSelectiveMode(select)
        self.setCoorndinateType(coordinate)
        self._atoms_ = []

        if not(filename is None):
            self.readPOSCAR(filename)

    def setComment(self, comment):
        self._comment_ = comment
        
    def getComment(self):
        return self._comment_

    def setLattice(self, vectors, latticeConstant = 1.0):
        """ set lattice
            vectors:         lattice vectors {vector array}
            latticeConstant:                 {number}
        """
        v1 = vectors[0]
        v2 = vectors[1]
        v3 = vectors[2]
        self._lattice_ = Lattice(v1, v2, v3, latticeConstant)

    def getLattice(self):
        return self._lattice_

    def addAtom(self, atom):
#        a = Atom(atom._symbol_,
#                 atom._xCoordinate_, atom._yCoordinate_, atom._zCoordinate_,
#                 atom._xDynamic_,    atom._yDynamic_,    atom._zDynamic_,
#                 atom._xDisplace_,   atom._yDisplace_,   atom._zDisplace_)
        a = Atom()
        a.copyAtom(atom)
        self._atoms_.append(a)

    def delAtom(self):
        pass

    def setAtomElement(self, index, elementSymbol):
#        self._atoms_[index].setElement(element)
        self._atoms_[index].setSymbol(elementSymbol)

    def setAtomCoordinate(self, index, xCoordinate, yCoordinate, zCoordinate):
        self._atoms_[index].setCoordinate(xCoordinate, yCoordinate, zCoordinate)

    def setAtomDynamic(self, index, xDynamic, yDynamic, zDynamic):
        self._atoms_[index].setDynamic(xDynamic, yDynamic, zDynamic)

    def setAtom(self, index, atom):
        self._atoms_[index] = atom

    def listAtom(self):
        for a in self._atoms_:
            a.showAtom()
    
    def _checkElements_(self):
#        tmpElementType = self._atoms_[0].getElement()
        tmpElementType = self._atoms_[0].getSymbol()
        tmpElementNumber = 0
        elements = []
        for a in self._atoms_:
#            e = a.getElement()
            e = a.getSymbol()
            if tmpElementType == e:
                tmpElementNumber += 1
            else:
                elements.append({tmpElementType: tmpElementNumber})
                tmpElementType = e
                tmpElementNumber = 1
        elements.append({tmpElementType: tmpElementNumber})
        return elements

#    def setElementsType(self, elements = ['A', 'B', 'C' ,'D']):
    def setElementsType(self, elements):
        es = self._checkElements_()
        n = 0
        for i in range(len(es) ):
            for j in range(es[i].values()[0] ):
#                self._atoms_[n].setElement(elements[i])
                self._atoms_[n].setSymbol(elements[i])
                n += 1

    def setSelectiveMode(self, mode):
        self._selectiveMode_ = mode
        
    def getSelectiveMode(self):
        return self._selectiveMode_
    
    def setCoorndinateType(self, coordinate):
        self._coorndinateType_ = coordinate
        
    def getCoorndinateType(self):
        return self._coorndinateType_

    def readPOSCAR(self, filename):
        f = open(filename)
        self.setComment(f.readline().rstrip())
        # setup lattice constant
        latticeConstant = float(f.readline().split()[0])
        vectors = []

        # setup lattice vector
        for i in range(3):
#            l = f.readline().rstrip()
            l = f.readline().split()
            tmpVec = Vector(float(l[0]), float(l[1]), float(l[2]) )
            vectors.append(tmpVec)
#            vectors.append((float(l.split()[0]), float(l.split()[1]), float(l.split()[2]) ) )
        self.setLattice(vectors, latticeConstant)

        # setup number of element type
        tmpAtomTypes = []
        tmpAtomNumbers = []
        l = f.readline().split()

        # for VASP5 POSCAR
        if l[0].isalpha():
            tmpAtomTypes = l
            l = f.readline().split()

        totalAtomNumber = 0
        for n in l:
            tmpAtomNumbers.append(int(n) )

        # setup selective mode and coordinate type
        l = f.readline().rstrip()
        if l.split()[0][0].upper() == "S":
            # setup coordinate type
            self._selectiveMode_ = l
            l = f.readline().rstrip()
            # cartesian coordinates
            if l.split()[0][0].upper() == 'C':
                self._coorndinateType_ = l
            # cartesian coordinates
            elif l.split()[0][0].upper() == 'K':
                self._coorndinateType_ = l
            # direct/fractional coordinates)
            elif l.split()[0][0].upper() == 'D':
                self._coorndinateType_ = l
            else:
                pass

        # setup atom coordinate
        for i in range(len(tmpAtomNumbers) ):
            if len(tmpAtomTypes) is not 0:
                symbol = tmpAtomTypes[i]
            else:
                symbol = str(i)
                pass
            for j in range(tmpAtomNumbers[i]):
                l = f.readline().split()
                if len(l) == 3:
                    a = Atom(symbol ,float(l[0]) ,float(l[1]), float(l[2]) )
                    checkElementByPeriodicTable(a, 'symbol')
                elif len(l) == 6:
                    a = Atom(symbol ,float(l[0]) ,float(l[1]), float(l[2]), l[3], l[4], l[5])
                    checkElementByPeriodicTable(a, 'symbol')
                else:
                    print "error format"
                self.addAtom(a)

        f.close()

    def writePOSCAR(self, filename = None):
        comment = ''
        numberOfAtom = ''
        atomOfType = ''

        for es in self._checkElements_():
            atomOfType += ' ' + es.keys()[0]
            numberOfAtom += ' ' + str(es.values()[0])

        l = self._lattice_.getVectors()
        v11, v12, v13 = l[0].getBasis()
        v21, v22, v23 = l[1].getBasis()
        v31, v32, v33 = l[2].getBasis()
        format1 = '''%s
%.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
'''
        output1 = format1 % (self._comment_ + comment, self._lattice_.getConstant(),
                             v11, v12, v13,
                             v21, v22, v23,
                             v31, v32, v33)


        if atomOfType.split()[0].istitle():
            output2 = atomOfType + "\n"
            output2 += numberOfAtom + "\n"
        else:
            output2 = numberOfAtom + "\n"
#        output2 = numberOfAtom + "\n"
#        output += numberOfAtom + "\n"

        output3 = ""
        if self._selectiveMode_ is None:
            output3 = "%s\n" % (self._coorndinateType_)
            for a in self._atoms_:
                output3 += "%+14.10f %+14.10f %+14.10f\n" % a.getCoordinate()
        else:
            output3 = "%s\n%s\n" % (self._selectiveMode_, self._coorndinateType_)
            for a in self._atoms_:
                output3 += "%+14.10f %+14.10f %+14.10f" % a.getCoordinate() + "  %s %s %s\n" % a.getDynamic()
        
        # setup output
        if filename == None:
            print output1 + output2 + output3
        else:
            f = open(filename, "wb")
            f.write(output1 + output2 + output3)
            f.close()
#        self.listAtom()

    def lineScan(self, distance, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ line scan
            distance:              {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
        ref_atm = self._atoms_[ref_indexes[0]]
        mot_atm = self._atoms_[mot_indexes[0]]
        x1, y1, z1 = ref_atm.getCoordinate()
        x2, y2, z2 = mot_atm.getCoordinate()
        vec = Vector(x2-x1, y2-y1, z2-z1)
        vec = vec.normalized()
    
        poscars = []
        tmpIndexes = mot_indexes + grp_indexes
    
        for i in xrange(nstep + 1):
            v = vec*i*distance
#            tmpPOSCAR = copy.deepcopy(poscar)
            tmpPOSCAR = POSCAR(comment = self._comment_, select = self._selectiveMode_, coordinate = self._coorndinateType_)
            tmpPOSCAR._lattice_ = copy.deepcopy(self._lattice_)
            tmpPOSCAR._atoms_ = copy.deepcopy(self._atoms_)
            x, y, z = v.getBasis()
    #        tmpPOSCAR.setAtomCoordinate(mot_indexes[0], x2+x, y2 + y, z2 + z)
    
            for j in tmpIndexes:
                tmpX, tmpY, tmpZ = self._atoms_[j].getCoordinate()
                tmpPOSCAR.setAtomCoordinate(j, tmpX + x, tmpY + y, tmpZ + z)
    #            tmpX, tmpY, tmpZ = poscar._atoms_[j-1].getCoordinate()
    #            tmpPOSCAR.setAtomCoordinate(j-1, tmpX+x, tmpY+y, tmpZ+z)
            poscars.append(tmpPOSCAR)
        return poscars

    def angleScan(self, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ angle scan
            angle:                 {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
    # ref atm, fix/basic atm, mot atm
        ref_atm = self._atoms_[ref_indexes[0] ]
        bas_atm = self._atoms_[ref_indexes[1] ]
        mot_atm = self._atoms_[mot_indexes[0] ]
    
        x1, y1, z1 = ref_atm.getCoordinate()
        x2, y2, z2 = bas_atm.getCoordinate()
        x3, y3, z3 = mot_atm.getCoordinate()
        vec1 = Vector(x1 - x2, y1 - y2, z1 - z2)
        vec2 = Vector(x3 - x2, y3 - y2, z3 - z2)
        normal_vector = vec1.cross(vec2)
    
        poscars = []
        tmpIndexes = mot_indexes + grp_indexes
    
        for i in xrange(nstep):
            a = i * angle
#            tmpPOSCAR = copy.deepcopy(poscar)
            tmpPOSCAR = POSCAR(comment = self._comment_, select = self._selectiveMode_, coordinate = self._coorndinateType_)
            tmpPOSCAR._lattice_ = copy.deepcopy(self._lattice_)
            tmpPOSCAR._atoms_ = copy.deepcopy(self._atoms_)
    
    #        for j in grp_indexes:
            for j in tmpIndexes:
#                tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
                tmpX, tmpY, tmpZ = self._atoms_[j].getCoordinate()
                tmpV = Vector(tmpX - x2, tmpY - y2, tmpZ - z2)
                tmpV = tmpV.rotate(normal_vector, a)
                tmpX, tmpY, tmpZ = tmpV.getBasis()
                tmpPOSCAR.setAtomCoordinate(j, tmpX + x2, tmpY + y2, tmpZ + z2)
            poscars.append(tmpPOSCAR)
        return poscars

    def dihedralScan(self, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ dihedral scan
            angle:                 {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
    # ref atm, fix/basic atm, mot atm
        ref1_atm = self._atoms_[ref_indexes[0] ]
        ref2_atm = self._atoms_[ref_indexes[1] ]
        ref3_atm = self._atoms_[ref_indexes[2] ]
        mot_atm  = self._atoms_[mot_indexes[0] ]
    
        x1, y1, z1 = ref1_atm.getCoordinate()
        x2, y2, z2 = ref2_atm.getCoordinate()
        x3, y3, z3 = ref3_atm.getCoordinate()
        x4, y4, z4 = mot_atm.getCoordinate()
    
        vec1 = Vector(x1 - x2, y1 - y2, z1 - z2)
        vec2 = Vector(x3 - x2, y3 - y2, z3 - z2)
        vec3 = Vector(x2 - x3, y2 - y3, z2 - z3)
        vec4 = Vector(x4 - x3, y4 - y3, z4 - z3)
    
        normal_vec1 = vec1.cross(vec2)
        normal_vec2 = vec3.cross(vec4)
        normal_vec3 = normal_vec1.cross(normal_vec2)
    
        poscars = []
        tmpIndexes = mot_indexes + grp_indexes
    
        for i in xrange(nstep):
            a = i * angle
#            tmpPOSCAR = copy.deepcopy(poscar)
            tmpPOSCAR = POSCAR(comment = self._comment_, select = self._selectiveMode_, coordinate = self._coorndinateType_)
            tmpPOSCAR._lattice_ = copy.deepcopy(self._lattice_)
            tmpPOSCAR._atoms_ = copy.deepcopy(self._atoms_)
    
            for j in tmpIndexes:
#                tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
                tmpX, tmpY, tmpZ = self._atoms_[j].getCoordinate()
                tmpV = Vector(tmpX - x3, tmpY - y3, tmpZ - z3)
                tmpV = tmpV.rotate(normal_vec3, a)
                tmpX, tmpY, tmpZ = tmpV.getBasis()
                tmpPOSCAR.setAtomCoordinate(j, tmpX + x3, tmpY + y3, tmpZ + z3)
            poscars.append(tmpPOSCAR)
        return poscars

    def addPOSCARcoordindate(self, atoms):
        pass

    def directToCartesian(self):
        """ direct coordinate convert to cartesian coordinate
            ref: http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_cartesian_coordinates
        """
        if self.getCoorndinateType()[0].upper() == 'D':
            self.setCoorndinateType('Cartesian')
            l = self.getLattice().getVectors()
            lc = self.getLattice().getConstant()
            v1 = l[0]
            v2 = l[1]
            v3 = l[2]
            alpha = v2.getAngle(v3)
            beta  = v1.getAngle(v3)
            gamma = v1.getAngle(v2)
            v1Len = v1.getLength()
            v2Len = v2.getLength()
            v3Len = v3.getLength()
            cosAlpha = math.cos(math.radians(alpha) )
            cosBeta  = math.cos(math.radians(beta) )
            cosGamma = math.cos(math.radians(gamma) )
            sinAlpha = math.sin(math.radians(alpha) )
            sinBeta  = math.sin(math.radians(beta) )
            sinGamma = math.sin(math.radians(gamma) )
            volume = math.sqrt(1 - cosAlpha*cosAlpha - cosBeta*cosBeta - cosGamma*cosGamma + 2*cosAlpha*cosBeta*cosGamma)

            for a in self._atoms_:
                x, y, z = a.getCoordinate()
                tmpX = v1Len*x + v2Len*cosGamma*y + v3Len*cosBeta*z
                tmpY = v2Len*sinGamma*y + v3Len*(cosAlpha - cosBeta*cosGamma)/sinGamma*z
                tmpZ = v3Len*volume/sinGamma*z
                a.setCoordinate(tmpX, tmpY, tmpZ)

    def cartesianToDirect(self):
        """ cartesian coordinate convert to direct coordinate
            ref: http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_from_cartesian_coordinates
        """
        if self.getCoorndinateType()[0].upper() == 'C':
            self.setCoorndinateType('Direct')
            l = self.getLattice().getVectors()
            lc = self.getLattice().getConstant()
            v1 = l[0]
            v2 = l[1]
            v3 = l[2]
            alpha = v2.getAngle(v3)
            beta  = v1.getAngle(v3)
            gamma = v1.getAngle(v2)
            v1Len = v1.getLength()
            v2Len = v2.getLength()
            v3Len = v3.getLength()
            cosAlpha = math.cos(math.radians(alpha) )
            cosBeta  = math.cos(math.radians(beta) )
            cosGamma = math.cos(math.radians(gamma) )
            sinAlpha = math.sin(math.radians(alpha) )
            sinBeta  = math.sin(math.radians(beta) )
            sinGamma = math.sin(math.radians(gamma) )
            volume = math.sqrt(1 - cosAlpha*cosAlpha - cosBeta*cosBeta - cosGamma*cosGamma + 2*cosAlpha*cosBeta*cosGamma)

            for a in self._atoms_:
                x, y, z = a.getCoordinate()
                tmpX = 1/v1Len*x - cosGamma/v1Len/sinGamma*y + (cosAlpha*cosGamma - cosBeta)/v1Len/volume/sinGamma*z
                tmpY = 1/v2Len/sinGamma*y + (cosBeta*cosGamma - cosAlpha)/v2Len/volume/sinGamma*z
                tmpZ = sinGamma/v3Len/volume*z
                a.setCoordinate(tmpX, tmpY, tmpZ)
#                a.setCoordinate(tmpX/lc, tmpY/lc, tmpZ/lc)
            


class OUTCAR:
    def __init__(self, filename = 'OUTCAR'):
        self._elements_ = []
        self._dynamicMatrixes_ = []
        self.readOUTCAR(filename)

    def readOUTCAR(self, filename):
        f = open(filename)
        l = ' '
        reSpace = re.compile('^\s+?$')
#        reSpace = re.compile('^\s+?$')
        rePOTCAR = re.compile('^\s+?POTCAR:\s+?(\w+)\s+?(\w+)\s+?(\w+)')
#        rePosition = re.compile(' position of ions in cartesian coordinates  (Angst):')
        rePosition = re.compile('^\s+?position of ions in cartesian coordinates')
#        rePosition2 = re.compile('^\s+?(\w+)\s+?(\w+)\s+?(\w+)')
        reDynMat = re.compile('Eigenvectors and eigenvalues of the dynamical matrix')
        reFinite = re.compile('Finite differences POTIM')

        totalAtomNumber = 0
        while l:
            l = f.readline()
            # Get potcar
            if rePOTCAR.search(l):
                r = rePOTCAR.match(l)
                e1, e2, e3 = r.groups()
                element = {'potential': e1, 'element': e2, 'date': e3}
                self._elements_.append(element)

            # Get atom position
            if rePosition.match(l):
                l = f.readline()
                while not reSpace.search(l):
                    l = f.readline()
                    totalAtomNumber += 1
                    

            # Get dynamical matrix
            if reDynMat.search(l):
                l = f.readline()
#                while not reDynMat.search(l):
                while not (reDynMat.search(l) or reFinite.search(l) ):
                    tmpArray = l.split()
                    # Get image freq
                    if len(tmpArray) == 10:
                        freq = {"THz": float(tmpArray[2]) * -1.0,
                                "2PiTHz": float(tmpArray[4]) * 1.0,
                                "cm-1": float(tmpArray[6]) * -1.0,
                                "meV": float(tmpArray[8]) * -1.0}
                        atoms = []
                        l = f.readline()
                        tmpArray = l.split()

                        while len(tmpArray) == 6:
                            if tmpArray[0] == 'X':
                                l = f.readline()
                                tmpArray = l.split()
                            else:
#                                tmpAtm = Atom(element = 0,
                                tmpAtm = Atom(elementSymbol = 'X',
                                              xCoordinate = float(tmpArray[0]), yCoordinate = float(tmpArray[1]), zCoordinate = float(tmpArray[2]),
                                              xDisplace = float(tmpArray[3]), yDisplace = float(tmpArray[4]), zDisplace = float(tmpArray[3]) )
                                l = f.readline()
                                tmpArray = l.split()
                                atoms.append(tmpAtm)
                        self._dynamicMatrixes_.append({"freq": freq, "atoms": atoms})

                    # Get real freq
                    elif len(tmpArray) == 11:
                        freq = {"THz": float(tmpArray[3]),
                                "2PiTHz": float(tmpArray[5]),
                                "cm-1": float(tmpArray[7]),
                                "meV": float(tmpArray[9])}
                        atoms = []
                        l = f.readline()
                        tmpArray = l.split()

                        while len(tmpArray) == 6:
                            if tmpArray[0] == 'X':
                                l = f.readline()
                                tmpArray = l.split()
                            else:
#                                tmpAtm = Atom(element = 0,
                                tmpAtm = Atom(elementSymbol = 'X',
                                              xCoordinate = float(tmpArray[0]), yCoordinate = float(tmpArray[1]), zCoordinate = float(tmpArray[2]),
                                              xDisplace = float(tmpArray[3]), yDisplace = float(tmpArray[4]), zDisplace = float(tmpArray[3]) )
                                l = f.readline()
                                tmpArray = l.split()
                                atoms.append(tmpAtm)
#                            print l.rstrip()
                        self._dynamicMatrixes_.append({"freq": freq, "atoms": atoms})
 
                    l = f.readline()
                
        f.close()

    def writeLog(self, filename = None):
        out0orientation = """                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic     Atomic              Coordinates (Angstroms)
 Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------\n"""
        out0coordinate = '     %3d       %3d        %3d       %10.6f   %10.6f %10.6f\n'
        out0dash = '---------------------------------------------------------------------\n'
        out0frequency = """ Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:\n"""
        out1number = '                   %3d\n'
        out2number = '                   %3d                    %3d\n'
        out3number = '                   %3d                    %3d                         %3d\n'
        out1frequency = ' Frequencies --   %10.4f\n'
        out2frequency = ' Frequencies --   %10.4f             %10.4f\n'
        out3frequency = ' Frequencies --   %10.4f             %10.4f             %10.4f\n'
        out1title = ' Atom AN      X      Y      Z\n'
        out2title = ' Atom AN      X      Y      Z        X      Y      Z\n'
        out3title = ' Atom AN      X      Y      Z        X      Y      Z        X      Y      Z\n'
        out1col = ' %3d %3d   %6.2f %6.2f %6.2f\n'
        out2col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f\n'
        out3col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f\n'
        sentances = [{'orientation': out0orientation, 'coordinates': out0coordinate, 'dash': out0dash, 'frequency': out0frequency},
                     {'number': out1number, 'frequency': out1frequency, 'title': out1title, 'col': out1col},
                     {'number': out2number, 'frequency': out2frequency, 'title': out2title, 'col': out2col},
                     {'number': out3number, 'frequency': out3frequency, 'title': out3title, 'col': out3col},]

        numberOfFreq = len(self._dynamicMatrixes_)
        quotient = numberOfFreq / 3
        remainder = numberOfFreq % 3
        numberOfAtoms = len(self._dynamicMatrixes_[0]['atoms'])
#        print numberOfFreq, quotient, remainder
        out = sentances[0]['orientation']

        for i in range(numberOfAtoms):
            x, y, z = self._dynamicMatrixes_[0]['atoms'][i].getCoordinate()
            out += sentances[0]['coordinates'] %(i+1, 1, 0, x, y, z)
        out += sentances[0]['dash']
 
        out += sentances[0]['frequency']
        for i in range(quotient):
            out += sentances[3]['number'] %(3*i+1, 3*i+2, 3*i+3)
            out += sentances[3]['frequency'] %(self._dynamicMatrixes_[3*i-3]['freq']['cm-1'], self._dynamicMatrixes_[3*i-2]['freq']['cm-1'], self._dynamicMatrixes_[3*i-1]['freq']['cm-1'])
            out += sentances[3]['title']
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[3*i-3]['atoms'][j].getDisplace()
                a21, a22, a23 = self._dynamicMatrixes_[3*i-2]['atoms'][j].getDisplace()
                a31, a32, a33 = self._dynamicMatrixes_[3*i-1]['atoms'][j].getDisplace()
                out += sentances[3]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23, a31, a32, a33)

        if remainder == 1:
            out += sentances[1]['number'] %(numberOfFreq)
            out += sentances[1]['frequency'] %(self._dynamicMatrixes_[numberOfFreq-1]['freq']['cm-1'])
            out += sentances[1]['title']
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[numberOfFreq-1]['atoms'][j].getDisplace()
                out += sentances[1]['col'] %(j+1, 1, a11, a12, a13)
        elif remainder == 2:
            out += sentances[2]['number'] %(numberOfFreq-1, numberOfFreq)
            out += sentances[2]['frequency'] %(self._dynamicMatrixes_[numberOfFreq-2]['freq']['cm-1'], self._dynamicMatrixes_[numberOfFreq-2]['freq']['cm-1'])
            out = out + sentances[2]['title']
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[numberOfFreq-2]['atoms'][j].getDisplace()
                a21, a22, a23 = self._dynamicMatrixes_[numberOfFreq-1]['atoms'][j].getDisplace()
                out += sentances[2]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23)
#        out += '  -------------------------------------------------------------------'
        # setup output
        if filename == None:
            print out
        else:
            f = open(filename, "wb")
            f.write(out)
            f.close()
#    


if __name__ == "__main__":
    import sys
    import os

#    o = OUTCAR('OUTCAR')
#    o.writeLog()

#    p = POSCAR('POSCAR5C')
    p = POSCAR('c1')
    p.cartesianToDirect()
    p.writePOSCAR('d1')

    p = POSCAR('d1')
    p.directToCartesian()
    p.writePOSCAR('c2')

    pass
