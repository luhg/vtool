#!/usr/bin/env python
#
### modified date: 2013/10/23
#

from atm import *
#import copy
import math
import os

#def calcVector(coordinate1, coordinate2):
#    x = coordinate1[0] - coordinate2[0]
#    y = coordinate1[1] - coordinate2[1]
#    z = coordinate1[2] - coordinate2[2]
#    return (x, y, z)

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
        self._length_ = self._x_*self._x_ + self._y_*self._y_ + self._z_*self._z_
        self._length_ = math.sqrt(self._length_)
        return self._length_

    def dot(self, vector):
        pass

    def cross(self, vector):
        pass

    def __str__(self):
        return "%f, %f, %f" %(self._x_, self._y_, self._z_)

    def __add__(self, vector):
        return Vector(self._x_+vector._x_, self._y_+vector._y_, self._z_+vector._z_)

    def __sub__(self, vector):
        return Vector(self._x_-vector._x_, self._y_-vector._y_, self._z_-vector._z_)

    def __mul__(self, n):
        return Vector(n*self._x_, n*self._y_, n*self._z_)


def lineScan(file, distance, nstep, ref_index, mot_index, grp_indexes):
    poscar = POSCAR(file)
    ref_atm = poscar._atoms_[ref_index]
    mot_atm = poscar._atoms_[mot_index]
    x1, y1, z1 = ref_atm.getCoordinate()
    x2, y2, z2 = mot_atm.getCoordinate()
    vec = Vector(x2-x1, y2-y1, z2-z1)
    vec = vec.normalized()

    poscars = []

    for i in xrange(nstep + 1):
        v = vec*i*m_distance
        x, y, z = v.getBasis()
        tmpPOSCAR = POSCAR(file)
#        tmpPOSCAR = copy.copy(poscar)
        tmpPOSCAR.setAtomCoordinate(mot_index, x2+x, y2 + y, z2 + z)

        for j in grp_indexes:
            tmpX, tmpY, tmpZ = poscar._atoms_[j-1].getCoordinate()
            tmpPOSCAR.setAtomCoordinate(j-1, tmpX+x, tmpY+y, tmpZ+z)

        print i, x, y, z
        poscars.append(tmpPOSCAR)
    return poscars

def makeScanJob(poscars):
    i = 0
    for p in poscars:
        dir = "0%d" %(i) if i < 10 else "%2d" %(i)
        os.mkdir(dir)
        print os.path.join(dir, 'POSCAR')
        p.writePOSCAR(os.path.join(dir, 'POSCAR'))
        i = i + 1
    pass


if __name__ == "__main__":
#    import sys
#    import os

#    file = sys.argv[1]
    file = "iPOSCAR"
#    p = POSCAR(file)
#    g = GJF()
#    poscar2gjf(p, g)

#    line split
#   r: reference atom
#   m: motion atom
#   m_distance: motion distance
#   nstep: number of step
#   input
    ref_index = int('1') - 1
    mot_index = int('2') - 1
    m_distance = -0.3 
    nstep = 3
    grp_indexes = [3,4,6] 

    poscars = lineScan(file, m_distance, nstep, ref_index, mot_index, grp_indexes)
    makeScanJob(poscars)


#    angle split
#    dihedral split
    pass
