#!/usr/bin/env python
#
### modified date: 2013/10/21
#

from operator import itemgetter, attrgetter
import copy

class Atom:
    def __init__(self, element,
                 xCoordinate, yCoordinate, zCoordinate,
                 xDynamic = 'T', yDynamic = 'T', zDynamic = 'T'):
        self.setElement(element)
        self.setCoordinate(xCoordinate, yCoordinate, zCoordinate)
        self.setDynamic(xDynamic, yDynamic, zDynamic)

    def setElement(self, element):
        self._element_ = element

    def getElement(self):
        return self._element_

    def setCoordinate(self, xCoordinate, yCoordinate, zCoordinate):
        self._xCoordinate_ = xCoordinate
        self._yCoordinate_ = yCoordinate
        self._zCoordinate_ = zCoordinate

    def getCoordinate(self):
        return (self._xCoordinate_, self._yCoordinate_, self._zCoordinate_)

    def setDynamic(self, xDynamic, yDynamic, zDynamic):
        self._xDynamic_ = xDynamic
        self._yDynamic_ = yDynamic
        self._zDynamic_ = zDynamic

    def getDynamic(self):
        return (self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def __repr__(self):  
        return repr((self._element_, self._xCoordinate_, self._yCoordinate_, self._zCoordinate_) )

    def showAtom(self):
        print "%3s, %+14.10f, %+14.10f, %+14.10f, %4s, %4s, %4s" %(self._element_, self._xCoordinate_, self._yCoordinate_, self._zCoordinate_, self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def addCoordinate(self, x, y, z):
        X = self._xCoordinate_ + x
        Y = self._yCoordinate_ + y
        Z = self._zCoordinate_ + z
        return Atom(self._element_, X, Y, Z, self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def subCoordinate(self, x, y, z):
        X = self._xCoordinate_ - x
        Y = self._yCoordinate_ - y
        Z = self._zCoordinate_ - z
        return Atom(self._element_, X, Y, Z, self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def mulCoordinate(self, f):
        X = self._xCoordinate_ * f
        Y = self._yCoordinate_ * f
        Z = self._zCoordinate_ * f
        return Atom(self._element_, X, Y, Z, self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def divCoordinate(self, f):
        X = self._xCoordinate_ / f
        Y = self._yCoordinate_ / f
        Z = self._zCoordinate_ / f
        return Atom(self._element_, X, Y, Z, self._xDynamic_, self._yDynamic_, self._zDynamic_)

#    def __add__(self, atom):
#        x, y, z = atom.getCoordinate()
#        X = self._xCoordinate_ + x
#        Y = self._yCoordinate_ + y
#        Z = self._zCoordinate_ + z
#        newAtom = Atom(self._element_, X, Y, Z, self._xDynamic_, self._yDynamic_, self._zDynamic_)
#        return newAtom

class Lattice():
    def __init__(self, v1, v2, v3, constant = 1.0):
        self.setVectors(v1, v2, v3)
        self.setConstant(constant)

    def setConstant(self, const):
        self._constant_ = const

    def getConstant(self):
        return self._constant_

    def setVectors(self, v1, v2, v3):
        self._Vectors_ = [v1, v2, v3]

    def getVectors(self):
        return self._Vectors_


class GJF:
    def __init__(self, filename = None, comment = None):
        self._specs_ = []
#        self._option_ = ''
        self.setOption()
        self.setComment()
        self.setCharge()
        self.setSpin()
        self._atoms_ = []
        self._elements_ = []
        self._element_lists_ = []
        self._numbers_ = []
#        self._lattices_ =  []
        if not(filename is None):
            self.readGJF(filename)

    def addAtom(self, atom):
        self._atoms_.append(atom)

    def getAtoms(self):
        return self._atoms_

    def getLattice(self):
        return self._lattice_

    def setLattice(self, vectors, latticeConstant = 1.0):
        v1 = vectors[0]
        v2 = vectors[1]
        v3 = vectors[2]
        self._lattice_ = Lattice(v1, v2, v3, latticeConstant)

    def sortAtoms(self):
#        self._atoms_ = sorted(self._atoms_, key=lambda atom: atom._element_)
        self._atoms_ = sorted(self._atoms_, key=lambda atom: (atom._element_, atom._zCoordinate_) )

    def getOption(self):
        return self._option_

    def setOption(self, opt = '# opt freq hf/3-21g'):
        self._option_ = opt

    def setComment(self, comment = "This line is comment"):
        self._comment_ = comment

    def getComment(self):
        return self._comment_

    def getCharge(self):
        return self._charge_

    def setCharge(self, charge = 0):
        self._charge_ = charge

    def getSpin(self):
        return self._spin_

    def setSpin(self, spin = 1):
        self._spin_ = spin

    def readGJF(self, filename):
        f = open(filename, "r")
        i = 0
        for l in f.readlines():
            if l[0] == "%":
                self._specs_.append(l.rstrip() )
                continue
            elif l[0] == "#":
#                self._option_ = l.rstrip()
                self.setOption(l.rstrip() )
                continue 
            elif l.strip() == "":
                continue
            elif len(l.split()) == 4:
                 atom = Atom(l.split()[0], float(l.split()[1]), float(l.split()[2]), float(l.split()[3]) )
                 self.addAtom(atom)
        self._checkLattice_()
        self.sortAtoms()

    def writeGJF(self, filename = None):

        output = ''
        for l in self._specs_:
            output += l + "\n"
        output += self._option_ + "\n\n"
        output += self._comment_ + "\n\n"
        output += '%i %i\n' %(self._charge_, self._spin_)
        for a in self._atoms_:
            output += "%-2s" % a.getElement() + "       %13.8f    %13.8f    %13.8f\n" % a.getCoordinate()

        # setup output
        if filename == None:
            print output
        else:
            f = open(filename, "wb")
            f.write(output)
            f.close()

    def _checkLattice_(self):
        i = 0
        lattice_indexes = []
        vectors = []
        for atom in self._atoms_:
            if atom.getElement() == "Tv":
#                (x, y, z) = atom.getCoordinate()
#                vectors.append([x, y, z])
                vectors.append(atom.getCoordinate() )
                lattice_indexes.append(i)
            i += 1
        lattice_indexes.reverse()
        for l in lattice_indexes:
            self._atoms_.pop(l)
#        self._lattice_ = Lattice(vectors[0], vectors[1], vectors[2])
        self.setLattice(vectors)

class POSCAR:
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

def gjf2poscar(gjf, poscar):
    poscar.setLattice(gjf._lattice_.getVectors() )
    for a in gjf._atoms_:
        poscar.addAtom(a)

#def poscar2gjf(poscar, gjf, elements = ['A', 'B', 'C', 'D']):
def poscar2gjf(poscar, gjf, elements = None):
    if elements != None:
        poscar.setElementsType(elements)
    for a in poscar._atoms_:
        gjf.addAtom(a)
    for v in poscar.getLattice().getVectors():
        a = Atom('Tv', v[0], v[1], v[2])
        gjf.addAtom(a)

#def calcVector(coordinate1, coordinate2):
#    x = coordinate1[0] - coordinate2[0]
#    y = coordinate1[1] - coordinate2[1]
#    z = coordinate1[2] - coordinate2[2]
#    return (x, y, z)
#def calcVector(init_atom, final_atom):
#    ic = init_atom.getCoordinate()
#    fc = final_atom.getCoordinate()
#    x = fc[0] - ic[0]
#    y = fc[1] - ic[1]
#    z = fc[2] - ic[2]
#    return (x, y, z)


if __name__ == "__main__":
    import sys
    import os

#    file = sys.argv[1]
    file = "iPOSCAR"
    p = POSCAR(file)
#    g = GJF()
#    poscar2gjf(p, g)

#    line split
    r_index = int('2') - 1
    x, y, z = (4.0, 3.0, 2.0)

    old_atm = p._atoms_[r_index]
    new_atm = copy.copy(old_atm)
#    new_atm.setElement('a')
    new_atm.setCoordinate(x, y, z)

#    print old_atm.getCoordinate()
    new_atm = new_atm.addCoordinate(x, y, z)
    print old_atm
    print new_atm
#    print calcVector(old_atm, new_atm)

    print old_atm
    p.setAtomElement(0, 'b')
    p.writePOSCAR()
#    angle split
#    dihedral split
    pass
