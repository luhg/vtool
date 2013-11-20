#!/usr/bin/env python
#
### modified date: 2013/10/21
#

from atm import *

__author__ = ""
__date__ = "2013/11/21"

__version__ = "$Revision: 0.1$"


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
        self._atom_[index] = atom

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
#        self.writePOSCAR()

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
            vectors.append((float(l.split()[0]), float(l.split()[1]), float(l.split()[2]) ) )
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

        pass

    def writePOSCAR(self, filename = None):
        comment = ''
        numberOfAtom = ''
        for es in self._checkElements_():
            comment += ' ' + es.keys()[0]
            numberOfAtom += ' ' + str(es.values()[0])

        l = self._lattice_.getVectors()
        v1 = l[0]
        v2 = l[1]
        v3 = l[2]
        format1 = '''%s
%.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
'''
        output1 = format1 % (self._comment_ + comment, self._lattice_.getConstant(),
                             l[0][0], l[0][1], l[0][2],
                             l[1][0], l[1][1], l[1][2],
                             l[2][0], l[2][1], l[2][2])

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

    def addPOSCARcoordindate(atoms):
        pass

if __name__ == "__main__":
    import sys
    import os

    pass
