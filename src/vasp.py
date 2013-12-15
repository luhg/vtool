#!/usr/bin/env python

from atm import *
import re
import copy

__author__ = ""
__date__ = "2013/11/25"

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
        a = Atom(atom._element_, atom._xCoordinate_, atom._yCoordinate_, atom._zCoordinate_, atom._xDynamic_, atom._yDynamic_, atom._zDynamic_)
        self._atoms_.append(a)

    def delAtom(self):
        pass

    def setAtomElement(self, index, element):
        self._atoms_[index].setElement(element)

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
        tmpElementType = self._atoms_[0].getElement()
        tmpElementNumber = 0
        elements = []
        for a in self._atoms_:
            e = a.getElement()
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
                self._atoms_[n].setElement(elements[i])
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
        latticeConstant = float(f.readline().rstrip()[0])
        vectors = []

        # setup lattice vector
        for i in range(3):
            l = f.readline().rstrip()
            tmpVec = Vector(float(l.split()[0]), float(l.split()[1]), float(l.split()[2]) )
            vectors.append(tmpVec)
#            vectors.append((float(l.split()[0]), float(l.split()[1]), float(l.split()[2]) ) )
        self.setLattice(vectors, latticeConstant)

        # setup number of element type
        l = f.readline()
        tmpAtomNumbers = []
        totalAtomNumber = 0
        for n in l.split():
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
            for j in range(tmpAtomNumbers[i]):
                l = f.readline().split()
                if len(l) == 3:
                    a = Atom(str(i) ,float(l[0]) ,float(l[1]))
                elif len(l) == 6:
                    a = Atom(str(i) ,float(l[0]) ,float(l[1]) ,float(l[2]), l[3], l[4], l[5])
                else:
                    print "error format"
                self.addAtom(a)

        f.close()

    def writePOSCAR(self, filename = None):
        comment = ''
        numberOfAtom = ''
        for es in self._checkElements_():
            comment += ' ' + es.keys()[0]
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


        output2 = numberOfAtom + "\n"
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
                tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
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
                tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
                tmpV = Vector(tmpX - x3, tmpY - y3, tmpZ - z3)
                tmpV = tmpV.rotate(normal_vec3, a)
                tmpX, tmpY, tmpZ = tmpV.getBasis()
                tmpPOSCAR.setAtomCoordinate(j, tmpX + x3, tmpY + y3, tmpZ + z3)
            poscars.append(tmpPOSCAR)
        return poscars

    def addPOSCARcoordindate(atoms):
        pass

class OUTCAR:
    def __init__(self, filename = 'OUTCAR'):
        self._elements_ = []
        self._dynamicMatrixes_ = []
        self.readOUTCAT(filename)

    def readOUTCAT(self, filename):
        f = open(filename)
        l = ' '
        reSpace = re.compile('^\s+?$')
#        reSpace = re.compile('^\s+?$')
        rePOTCAR = re.compile('^\s+?POTCAR:\s+?(\w+)\s+?(\w+)\s+?(\w+)')
#        rePosition = re.compile(' position of ions in cartesian coordinates  (Angst):')
        rePosition = re.compile('^\s+?position of ions in cartesian coordinates')
#        rePosition2 = re.compile('^\s+?(\w+)\s+?(\w+)\s+?(\w+)')
        reDynMat = re.compile('Eigenvectors and eigenvalues of the dynamical matrix')

#        for l in f.readlines():
#            print l
        totalAtomNumber = 0
        while l:
            l = f.readline()
            # Get potcar
            if rePOTCAR.search(l):
                r = rePOTCAR.match(l)
                e1, e2, e3 = r.groups()
#                print l, e1, e2, e3
                element = {'potential': e1, 'element': e2, 'date': e3}
                self._elements_.append(element)

            # Get atom position
            if rePosition.match(l):
#                print l.rstrip()
                l = f.readline()
                while not reSpace.search(l):
                    l = f.readline()
#                    print l.rstrip(), totalAtomNumber
                    totalAtomNumber += 1
                    

            # Get dynamical matrix
            if reDynMat.search(l):
                l = f.readline()
                while not reDynMat.search(l):
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
                                tmpAtm = Atom(element = 0,
                                              xCoordinate = float(tmpArray[0]), yCoordinate = float(tmpArray[1]), zCoordinate = float(tmpArray[2]),
                                              xDisplace = float(tmpArray[3]), yDisplace = float(tmpArray[4]), zDisplace = float(tmpArray[3]) )
                                l = f.readline()
                                tmpArray = l.split()
                                atoms.append(tmpAtm)
#                            print l.rstrip()
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
                                tmpAtm = Atom(element = 0,
                                              xCoordinate = float(tmpArray[0]), yCoordinate = float(tmpArray[1]), zCoordinate = float(tmpArray[2]),
                                              xDisplace = float(tmpArray[3]), yDisplace = float(tmpArray[4]), zDisplace = float(tmpArray[3]) )
                                l = f.readline()
                                tmpArray = l.split()
                                atoms.append(tmpAtm)
#                            print l.rstrip()
                        self._dynamicMatrixes_.append({"freq": freq, "atoms": atoms})
 
                    l = f.readline()
                
        f.close()
#        print totalAtomNumber
        pass
    def writeOUTCAR(self):
#        numberOfFrequency = len(self._frequencies_)
        out1number = '                   %3d'
        out2number = '                   %3d                    %3d'
        out3number = '                   %3d                    %3d                         %3d'
        out1frequency = ' Frequencies --   %10.4f'
        out2frequency = ' Frequencies --   %10.4f             %10.4f'
        out3frequency = ' Frequencies --   %10.4f             %10.4f             %10.4f'
        out1title = ' Atom AN      X      Y      Z'
        out2title = ' Atom AN      X      Y      Z        X      Y      Z'
        out3title = ' Atom AN      X      Y      Z        X      Y      Z        X      Y      Z'
        out1col = ' %3d %3d   %6.2f %6.2f %6.2f'
        out2col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f'
        out3col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f'
        sentances = [{'number': out3number, 'frequency': out3frequency, 'title': out3title, 'col': out3col},
                     {'number': out1number, 'frequency': out1frequency, 'title': out1title, 'col': out1col},
                     {'number': out2number, 'frequency': out2frequency, 'title': out2title, 'col': out2col}]

        numberOfFreq = len(self._dynamicMatrixes_)
        quotient = numberOfFreq / 3
        remainder = numberOfFreq % 3
        numberOfAtoms = len(self._dynamicMatrixes_[0]['atoms'])
        print numberOfFreq, quotient, remainder
        for i in range(quotient):
            out1 = sentances[0]['number'] %(3*i+1, 3*i+2, 3*i+3)
            out2 = sentances[0]['frequency'] %(self._dynamicMatrixes_[3*i-3]['freq']['cm-1'], self._dynamicMatrixes_[3*i-2]['freq']['cm-1'], self._dynamicMatrixes_[3*i-1]['freq']['cm-1'])
            out3 = sentances[0]['title']
            print out1
            print out2
            print out3
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[3*i-3]['atoms'][j].getDisplace()
                a21, a22, a23 = self._dynamicMatrixes_[3*i-2]['atoms'][j].getDisplace()
                a31, a32, a33 = self._dynamicMatrixes_[3*i-1]['atoms'][j].getDisplace()
                out4 = sentances[0]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23, a31, a32, a33)
                print out4

        if remainder == 1:
            out1 = sentances[1]['number'] %(numberOfFreq)
            out2 = sentances[1]['frequency'] %(self._dynamicMatrixes_[numberOfFreq-1]['freq']['cm-1'])
            out3 = sentances[1]['title']
            print out1
            print out2
            print out3
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[numberOfFreq-1]['atoms'][j].getDisplace()
                out4 = sentances[1]['col'] %(j+1, 1, a11, a12, a13)
                print out4
        elif remainder == 2:
            out1 = sentances[2]['number'] %(numberOfFreq-1, numberOfFreq)
            out2 = sentances[2]['frequency'] %(self._dynamicMatrixes_[numberOfFreq-2]['freq']['cm-1'], self._dynamicMatrixes_[numberOfFreq-2]['freq']['cm-1'])
            out3 = sentances[2]['title']
            print out1
            print out2
            print out3
            for j in range(numberOfAtoms):
                a11, a12, a13 = self._dynamicMatrixes_[numberOfFreq-2]['atoms'][j].getDisplace()
                a21, a22, a23 = self._dynamicMatrixes_[numberOfFreq-1]['atoms'][j].getDisplace()
                out4 = sentances[2]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23)
                print out4
             

if __name__ == "__main__":
    import sys
    import os

#    o = OUTCAR('OUTCAR')
    o = OUTCAR('O2')
    o.writeOUTCAR()
#    p = POSCAR('POSCAR')
#    p.writePOSCAR()
    pass
