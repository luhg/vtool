#!/usr/bin/env python
#
### modified date: 2013/11/21
#

from operator import itemgetter, attrgetter
import math

__author__ = ""
__date__ = "2013/11/21"
__version__ = "$Revision: 0.1$"


class Vector:
    def __init__(self, x, y, z):
        self.setBasis(x, y, z)

    def setBasis(self, x, y, z):
        self._x_ = x
        self._y_ = y
        self._z_ = z

    def getBasis(self):
        return (self._x_, self._y_, self._z_)

    def normalized(self):
        length = self.getLength()
        if length == 0.0:
            print "Don't normalized"
            return Vector(self._x_, self._y_, self._z_)
        else:
            return Vector(self._x_ / length, self._y_ / length, self._z_ / length)

    def getLength(self):
        self._length_ = self._x_ * self._x_ + self._y_ * self._y_ + self._z_ * self._z_
        self._length_ = math.sqrt(self._length_)
        return self._length_

    def getAngle(self, vector):
        angle = math.acos(self.dot(vector) / self.getLength() / vector.getLength() )
        return math.degrees(angle)

    def dot(self, vector):
        return self._x_ * vector._x_ + self._y_ * vector._y_ + self._z_ * vector._z_

    def cross(self, vector):
        x = self._y_ * vector._z_ - self._z_ * vector._y_
        y = self._z_ * vector._x_ - self._x_ * vector._z_
        z = self._x_ * vector._y_ - self._y_ * vector._x_
        return Vector(x, y, z)

    def __str__(self):
        return "%f, %f, %f" %(self._x_, self._y_, self._z_)

    def __add__(self, vector):
        return Vector(self._x_ + vector._x_, self._y_ + vector._y_, self._z_ + vector._z_)

    def __sub__(self, vector):
        return Vector(self._x_ - vector._x_, self._y_ - vector._y_, self._z_ - vector._z_)

    def __mul__(self, n):
        return Vector(n * self._x_, n * self._y_, n * self._z_)

# Rotation operator
#      cosQ+Nx^2(1-cosQ)        NxNy(1-cosQ) - Nz*sinQ   NxNz(1-cosQ) + Ny*sinQ     x1   x2   x3
# R =  NyNz(1-cosQ) + Nz*sinQ   cosQ + Ny^2(1-cosQ)      NyNz(1-cosQ) - Nz*sinQ  =  y1   y2   y3
#      NzNx(1-cosQ) - Ny*sinQ   NzNy(1-cosQ) + NxsinQ    cosQ + Nz^2(1-csoQ)a       z1   z2   z3
#
# angle degree to radian
#    def rotate(axis_vector, angle):
    def rotate(self, axis_vector, angle):
        angle = math.radians(angle)
        x1 = math.cos(angle) + axis_vector._x_ * axis_vector._x_ * (1 -math.cos(angle) )
        x2 = axis_vector._x_ * axis_vector._y_ * (1 - math.cos(angle) ) - axis_vector._z_ * math.sin(angle)
        x3 = axis_vector._x_ * axis_vector._z_ * (1 - math.cos(angle) ) + axis_vector._y_ * math.sin(angle)

        y1 = axis_vector._y_ * axis_vector._z_ * (1- math.cos(angle) ) + axis_vector._z_ * math.sin(angle)
        y2 = math.cos(angle) + axis_vector._y_ * axis_vector._y_ * (1- math.cos(angle) )
        y3 = axis_vector._y_ * axis_vector._z_ * (1- math.cos(angle) ) - axis_vector._z_ * math.sin(angle)

        z1 = axis_vector._z_ * axis_vector._x_ * (1- math.cos(angle) ) - axis_vector._y_ * math.sin(angle)
        z2 = axis_vector._z_ * axis_vector._y_ * (1- math.cos(angle) ) + axis_vector._x_ * math.sin(angle)
        z3 = math.cos(angle) + axis_vector._z_ * axis_vector._z_ * (1- math.cos(angle) )

        x = x1 * self._x_ + x2 * self._y_ + x3 * self._z_
        y = y1 * self._x_ + y2 * self._y_ + y3 * self._z_
        z = z1 * self._x_ + z2 * self._y_ + z3 * self._z_
        return Vector(x, y, z)


class Lattice():
    def __init__(self, v1, v2, v3, constant = 1.0):
        """ set default argment
            v1:       lattice vector1  {Vector}
            v2:       lattice vector2  {Vector}
            v3:       lattice vector3  {Vector}
            constant: lattice constant {Number}
        """
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


class Element:
    def __init__(self, symbol = 'Z', name = 'Dummy', number = 0, mass = 0.0):
        """ set default
            symbol:   chemical symbol {String}  [Z]
            name:     name of element {String}  [Dummy]
            number:   atomic number   {Int}     [0]
            mass:     atomic mass     {Float}   [0.0]
        """
        self.setSymbol(symbol)
        self.setName(name)
        self.setAtomicNumber(number)
        self.setAtomicMass(mass)

    def setSymbol(self, s = 'Z'):
        self._symbol_ = s

    def setName(self, n = 'Dummy'):
        self._name_ = n

    def setAtomicNumber(self, an = 0):
        self._atomicNumber_ = an

    def setAtomicMass(self, am = 0.0):
        self._atomicMass_ = am

    def getSymbol(self):
        return self._symbol_

    def getName(self):
        return self._name_

    def getAtomicNumber(self):
        return self._atomicNumber_

    def getAtomicMass(self):
        return self._atomicMass_

    def copyElement(self, e):
        self.setSymbol(e._symbol_)
        self.setName(e._name_)
        self.setAtomicNumber(e._atomicNumber_)
        self.setAtomicMass(e._atomicMass_)


PERIODIC_TABLE_ElEMENTS = [Element(),
                           Element('H',  'Hydrogen', 1, 1.00794),
                           Element('He', 'Helium',   2, 4.002602),
                           Element('Li', 'Lithium',  3, 6.941),
                           ]


def checkElementByPeriodicTable(element, method = 'symbol'):
    if method == 'symbol':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.getSymbol() == element.getSymbol():
                element.copyElement(e)
                break
    elif method == 'name':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.getName() == element.getName():
                element.copyElement(e)
                break

class Atom:
    """ Atom basic info """
    def __init__(self, element,
                 xCoordinate, yCoordinate, zCoordinate,
                 xDynamic = 'T', yDynamic = 'T', zDynamic = 'T',
                 xDisplace = 0.0, yDisplace = 0.0, zDisplace = 0.0):
        """ set default argments
            elelment:    atomic element
            xCoordinate: x-axix coordinate
            yCoordinate: y-axix coordinate
            zCoordinate: z-axix coordinate
            xDynamic:    x-axix (T)ranslate/(F)reeze
            yDynamic:    y-axix (T)ranslate/(F)reeze
            zDynamic:    z-axix (T)ranslate/(F)reeze
            xDisplace:   x-axix displacement
            yDisplace:   y-axix displacement
            zDisplace:   z-axix displacement
        """
        self.setElement(element)
        self.setCoordinate(xCoordinate, yCoordinate, zCoordinate)
        self.setDynamic(xDynamic, yDynamic, zDynamic)
        self.setDispalce(xDisplace, yDisplace, zDisplace)

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

    def setDispalce(self, xDisplace, yDisplace, zDisplace):
        self._xDisplace_ = xDisplace
        self._yDisplace_ = yDisplace
        self._zDisplace_ = zDisplace

    def getDisplace(self):
        return (self._xDisplace_, self._yDisplace_, self._zDisplace_)

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


if __name__ == "__main__":
    import sys
    import os

    pass
