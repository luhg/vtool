#!/usr/bin/env python
#
### modified date: 2013/11/21
#

from operator import itemgetter, attrgetter
import math

__author__ = ""
__date__ = "2013/12/19"
__version__ = "$Revision: 0.2$"


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

    def rotate(self, axis_vector, angle):
        """ rotating vector around arbitrary axis
                 cosQ+Nx^2(1-cosQ)        NxNy(1-cosQ) - Nz*sinQ   NxNz(1-cosQ) + Ny*sinQ     x1   x2   x3
            R =  NyNz(1-cosQ) + Nz*sinQ   cosQ + Ny^2(1-cosQ)      NyNz(1-cosQ) - Nz*sinQ  =  y1   y2   y3
                 NzNx(1-cosQ) - Ny*sinQ   NzNy(1-cosQ) + NxsinQ    cosQ + Nz^2(1-csoQ)a       z1   z2   z3
            ref: http://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis-angle
        """
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
    def __init__(self, symbol = 'X', name = 'Dummy', number = 0, mass = 0.0):
        """ set default
            symbol:   chemical symbol {String}  [X]
            name:     name of element {String}  [Dummy]
            number:   atomic number   {Int}     [0]
            mass:     atomic mass     {Float}   [0.0]
        """
        self.setSymbol(symbol)
        self.setName(name)
        self.setAtomicNumber(number)
        self.setAtomicMass(mass)

    def setSymbol(self, s = 'X'):
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


class Atom(Element):
    """ Atom basic info """
    def __init__(self, elementSymbol = 'X',
                 xCoordinate = 0.0, yCoordinate = 0.0, zCoordinate = 0.0,
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
#        self.setElement(element)
        self.setSymbol(elementSymbol)
        self.setCoordinate(xCoordinate, yCoordinate, zCoordinate)
        self.setDynamic(xDynamic, yDynamic, zDynamic)
        self.setDispalce(xDisplace, yDisplace, zDisplace)

#    def setElement(self, element):
#        self._element_ = element

#    def getElement(self):
#        return self._element_

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
#        return repr((self._element_, self._xCoordinate_, self._yCoordinate_, self._zCoordinate_) )
        return repr((self._symbol_, self._xCoordinate_, self._yCoordinate_, self._zCoordinate_) )

    def showAtom(self):
        print "%3s, %+14.10f, %+14.10f, %+14.10f, %4s, %4s, %4s" %(self._element_, self._xCoordinate_, self._yCoordinate_, self._zCoordinate_, self._xDynamic_, self._yDynamic_, self._zDynamic_)

    def copyAtom(self, atom):
        self.copyElement(atom)
#        self.setSymbol(atom._symbol_)
        self.setCoordinate(atom._xCoordinate_, atom._yCoordinate_, atom._zCoordinate_)
        self.setDynamic(atom._xDynamic_, atom._yDynamic_, atom._zDynamic_)
        self.setDispalce(atom._xDisplace_, atom._yDisplace_, atom._zDisplace_)
#        self.set
        pass

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


""" default periodic table element:
    global constant: PERIODIC_TABLE_ElEMENTS
"""
PERIODIC_TABLE_ElEMENTS = [Element(),
                           Element('H',   'Hydrogen',      1,   1.00794),
                           Element('He',  'Helium',        2,   4.002602),
                           Element('Li',  'Lithium',       3,   6.941),
                           Element('Be',  'Beryllium',     4,   9.0121831),
                           Element('B',   'Boron',         5,   10.81),
                           Element('C',   'Carbon',        6,   12.011),
                           Element('N',   'Nitrogen',      7,   14.007),
                           Element('O',   'Oxygen',        8,   15.999),
                           Element('F',   'Fluorine',      9,   18.998403163),
                           Element('Ne',  'Neon',          10,  20.1797),
                           Element('Na',  'Sodium',        11,  22.98976928),
                           Element('Mg',  'Magnesium',     12,  24.305),
                           Element('Al',  'Aluminium',     13,  26.9815385),
                           Element('Si',  'Silicon',       14,  28.085),
                           Element('P',   'Phosphorus',    15,  30.973761998),
                           Element('S',   'Sulfur',        16,  32.066),
                           Element('Cl'   'Chlorine',      17,  35.45),
                           Element('Ar',  'Argon',         18,  39.948),
                           Element('K',   'Potassium',     19,  39.0983),
                           Element('Ca',  'Calcium',       20,  40.078),
                           Element('Sc',  'Scandium',      21,  44.955908),
                           Element('Ti',  'Titanium',      22,  47.867),
                           Element('V',   'Vanadium',      23,  50.9415),
                           Element('Cr',  'Chromium',      24,  51.9961),
                           Element('Mn',  'Manganese',     25,  54.938044),
                           Element('Fe',  'Iron',          26,  55.845),
                           Element('Co',  'Cobalt',        27,  58.933194),
                           Element('Ni',  'Nickel',        28,  58.6934),
                           Element('Cu',  'Copper',        29,  63.546),
                           Element('Zn',  'Zinc',          30,  65.38),
                           Element('Ga',  'Gallium',       31,  69.723),
                           Element('Ge',  'Germanium',     32,  72.630),
                           Element('As',  'Arsenic',       33,  74.921595),
                           Element('Se',  'Selenium',      34,  78.971),
                           Element('Br',  'Bromine',       35,  79.904),
                           Element('Kr',  'Krypton',       36,  83.798),
                           Element('Rb',  'Rubidium',      37,  85.4678),
                           Element('Sr',  'Strontium',     38,  87.62),
                           Element('Y',   'Yttrium',       39,  88.90584),
                           Element('Zr',  'Zirconium',     40,  91.224),
                           Element('Nb',  'Niobium',       41,  92.90637),
                           Element('Mo',  'Molybdenum',    42,  95.95),
                           Element('Tc',  'Technetium',    43,  98),
                           Element('Ru',  'Ruthenium',     44,  101.07),
                           Element('Rh',  'Rhodium',       45,  102.90550),
                           Element('Pd',  'Palladium',     46,  106.42),
                           Element('Ag',  'Silver',        47,  107.8682),
                           Element('Cd',  'Cadmium',       48,  112.414),
                           Element('In',  'Indium',        49,  114.818),
                           Element('Sn',  'Tin',           50,  118.710),
                           Element('Sb',  'Antimony',      51,  121.760),
                           Element('Te',  'Tellurium',     52,  127.60),
                           Element('I',   'Iodine',        53,  126.90447),
                           Element('Xe',  'Xenon',         54,  131.293),
                           Element('Cs',  'Caesium',       55,  132.90545196),
                           Element('Ba',  'Barium',        56,  137.327),
                           Element('La',  'Lanthanum',     57,  138.90547),
                           Element('Ce',  'Cerium',        58,  140.116),
                           Element('Pr',  'Praseodymium',  59,  140.90766),
                           Element('Nd',  'Neodymium',     60,  144.242),
                           Element('Pm',  'Promethium',    61,  145),
                           Element('Sm',  'Samarium',      62,  150.36),
                           Element('Eu',  'Europium',      63,  151.964),
                           Element('Gd',  'Gadolinium',    64,  157.25),
                           Element('Tb',  'Terbium',       65,  158.92535),
                           Element('Dy',  'Dysprosium',    66,  162.500),
                           Element('Ho',  'Holmium',       67,  164.93033),
                           Element('Er',  'Erbium',        68,  167.259),
                           Element('Tm',  'Thulium',       69,  168.93422),
                           Element('Yb',  'Ytterbium',     70,  173.054),
                           Element('Lu',  'Lutetium',      71,  174.9668),
                           Element('Hf',  'Hafnium',       72,  178.49),
                           Element('Ta',  'Tantalum',      73,  180.94788),
                           Element('W',   'Tungsten',      74,  183.84),
                           Element('Re',  'Rhenium',       75,  186.207),
                           Element('Os',  'Osmium',        76,  190.23),
                           Element('Ir',  'Iridium',       77,  192.217),
                           Element('Pt',  'Platinum',      78,  195.084),
                           Element('Au',  'Gold',          79,  196.966569),
                           Element('Hg',  'Mercury',       80,  200.592),
                           Element('Tl',  'Thallium',      81,  204.38),
                           Element('Pb',  'Lead',          82,  207.2),
                           Element('Bi',  'Bismuth',       83,  208.98040),
                           Element('Po',  'Polonium',      84,  209),
                           Element('At',  'Astatine',      85,  210),
                           Element('Rn',  'Radon',         86,  222),
                           Element('Fr',  'Francium',      87,  223),
                           Element('Ra',  'Radium',        88,  226),
                           Element('Ac',  'Actinium',      89,  227),
                           Element('Th',  'Thorium',       90,  232.0377),
                           Element('Pa',  'Protactinium',  91,  231.03588),
                           Element('U',   'Uranium',       92,  238.02891),
                           Element('Np',  'Neptunium',     93,  237),
                           Element('Pu',  'Plutonium',     94,  244),
                           Element('Am',  'Americium',     95,  243),
                           Element('Cm',  'Curium',        96,  247),
                           Element('Bk',  'Berkelium',     97,  247),
                           Element('Cf',  'Californium',   98,  251),
                           Element('Es',  'Einsteinium',   99,  252),
                           Element('Fm',  'Fermium',       100, 257),
                           Element('Md',  'Mendelevium',   101, 258),
                           Element('No',  'Nobelium',      102, 259),
                           Element('Lr',  'Lawrencium',    103, 262),
                           Element('Rf',  'Rutherfordium', 104, 267),
                           Element('Db',  'Dubnium',       105, 268),
                           Element('Sg',  'Seaborgium',    106, 269),
                           Element('Bh',  'Bohrium',       107, 270),
                           Element('Hs',  'Hassium',       108, 269),
                           Element('Mt',  'Meitnerium',    109, 278),
                           Element('Ds',  'Darmstadtium',  110, 281),
                           Element('Rg',  'Roentgenium',   111, 281),
                           Element('Cn',  'Copernicium',   112, 285),
                           Element('Uut', 'Ununtrium',     113, 286),
                           Element('Fl',  'Flerovium',     114, 289),
                           Element('Uup', 'Ununpentium',   115, 288),
                           Element('Lv',  'Livermorium',   116, 293),
                           Element('Uus', 'Ununseptium',   117, 294),
                           Element('Uuo', 'Ununoctium',    118, 294),
                           Element('Tv',  'Translation vectors', -2, 0), # Defined in Gaussian
                          ]

def checkElementByPeriodicTable(element, method = 'symbol'):
    """ setup element all properties in peridic table
    """
    if method == 'symbol':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.getSymbol() == element.getSymbol():
                element.copyElement(e)
                break
            element.setName('dummy')
            element.setAtomicNumber(0)
            element.setAtomicMass(0.0)
    elif method == 'name':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.getName() == element.getName():
                element.copyElement(e)
                break
            element.setSymbol('X')
            element.setAtomicNumber(0)
            element.setAtomicMass(0.0)
    elif method == 'number':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.getAtomicNumber() == element.getAtomicNumber():
                element.copyElement(e)
                break
            element.setSymbol('X')
            element.setName('dummy')
            element.setAtomicMass(0.0)


if __name__ == "__main__":
    import sys
    import os

    pass
