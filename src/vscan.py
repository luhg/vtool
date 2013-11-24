#!/usr/bin/env python
#
### modified date: 2013/10/23
#

from atm import *
import copy
import os, getopt
import sys


def lineScan(poscar, distance, nstep, ref_indexes, mot_indexes, grp_indexes):
    ref_atm = poscar._atoms_[ref_indexes[0]]
    mot_atm = poscar._atoms_[mot_indexes[0]]
    x1, y1, z1 = ref_atm.getCoordinate()
    x2, y2, z2 = mot_atm.getCoordinate()
    vec = Vector(x2-x1, y2-y1, z2-z1)
    vec = vec.normalized()

    poscars = []
    tmpIndexes = mot_indexes + grp_indexes

    for i in xrange(nstep + 1):
        v = vec*i*distance
        tmpPOSCAR = cpy.deepcopy(poscar)
#        x, y, z = v.getBasis()
#        tmpPOSCAR.setAtomCoordinate(mot_indexes[0], x2+x, y2 + y, z2 + z)

        for j in grp_indexes:
            tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
            tmpPOSCAR.setAtomCoordinate(j, tmpX+x, tmpY+y, tmpZ+z)
#            tmpX, tmpY, tmpZ = poscar._atoms_[j-1].getCoordinate()
#            tmpPOSCAR.setAtomCoordinate(j-1, tmpX+x, tmpY+y, tmpZ+z)

#        print i, x, y, z
        poscars.append(tmpPOSCAR)
    return poscars


def angleScan(poscar, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
# ref atm, fix/basic atm, mot atm
    ref_atm = poscar._atoms_[ref_indexes[0] ]
    bas_atm = poscar._atoms_[ref_indexes[1] ]
    mot_atm = poscar._atoms_[mot_indexes[0] ]

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
        tmpPOSCAR = copy.deepcopy(poscar)

#        for j in grp_indexes:
        for j in tmpIndexes:
            tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
            tmpV = Vector(tmpX - x2, tmpY - y2, tmpZ - z2)
            tmpV = tmpV.rotate(normal_vector, a)
            tmpX, tmpY, tmpZ = tmpV.getBasis()
            tmpPOSCAR.setAtomCoordinate(j, tmpX + x2, tmpY + y2, tmpZ + z2)
        poscars.append(tmpPOSCAR)
    return poscars


def dihedralScan(poscar, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
# ref atm, fix/basic atm, mot atm
    ref1_atm = poscar._atoms_[ref_indexes[0] ]
    ref2_atm = poscar._atoms_[ref_indexes[1] ]
    ref3_atm = poscar._atoms_[ref_indexes[2] ]
    mot_atm = poscar._atoms_[mot_indexes[0] ]

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
        tmpPOSCAR = copy.deepcopy(poscar)

        for j in tmpIndexes:
            tmpX, tmpY, tmpZ = poscar._atoms_[j].getCoordinate()
            tmpV = Vector(tmpX - x3, tmpY - y3, tmpZ - z3)
            tmpV = tmpV.rotate(normal_vec3, a)
            tmpX, tmpY, tmpZ = tmpV.getBasis()
            tmpPOSCAR.setAtomCoordinate(j, tmpX + x3, tmpY + y3, tmpZ + z3)
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


def main():
    def usage():
        print "Usage: vscan -i POSCAR -r 1     -m 2 -g 3,4,5 -n 4 -d 1.0     # line  scan"
        print "       vscan -i POSCAR -r 1,2   -m 3 -g 4,5,6 -n 4 -d 30.0    # angle scan"
        print "       vscan -i POSCAR -r 1,2,3 -m 4 -g 5,6,7 -n 4 -d 30.0    # dihedral scan"
        print " -h : help"
        print " -i : input file, POSCAR"
        print " -r : reference atoms"
        print " -m : motion atom"
        print " -g : atom with motion atom"
        print " -n : number of step"
        print " -d : displacement (distance/angle)"

    def checkFile(f):
        if os.path.exists(f):
            return os.path.abspath(f)
        print "File: " + f + " didn't exist."
        sys.exit(2)
        return False

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], "hi:r:m:g:n:d:")

    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    infile = None
    r_indexes = None
    m_indexes = None
    g_indexes = None
    nstep = None
    displacement = None
    for o, a in opt_list:
        if o in ('-h'):
            usage()
            sys.exit()
        elif o in ('-i'):
            infile = a
        elif o in ('-r'):
#            r_indexes = a
            r_indexes = [int(i)-1 for i in a.split(',')]
        elif o in ('-m'):
#            m_indexes = a
            m_indexes = [int(i)-1 for i in a.split(',')]
        elif o in ('-g'):
#            g_indexes = a
            g_indexes = [int(i)-1 for i in a.split(',')]
        elif o in ('-n'):
            nstep = int(a)
        elif o in ('-d'):
            displacement = float(a)

    if infile is None:
        print "Intput file: "
        infile = sys.stdin.readline().rstrip()
#        infile = checkFile(sys.stdin.readline().rstrip())
    poscar = POSCAR(infile)

    if r_indexes is None:
        print "reference atoms: "
        r_indexes = sys.stdin.readline().rstrip()
        r_indexes = [int(i)-1 for i in r_indexes.split(',')]

    if m_indexes is None:
        print "motion atom: "
        m_indexes = sys.stdin.readline().rstrip()
        m_indexes = [int(i)-1 for i in m_indexes.split(',')]

    if g_indexes is None:
#        print "group atoms: "
#        g_indexes = sys.stdin.readline().rstrip()
#        g_indexes = [int(i)-1 for i in g_indexes.split(',')]
        g_indexes = []

    print r_indexes, m_indexes, g_indexes

    if nstep is None:
        print "number of step: "
        nstep = sys.stdin.readline().rstrip()
        nstep = int(nstep)

    if displacement is None:
        print "displacement distance/angle: "
        displacement = sys.stdin.readline().rstrip()
        displacement = float(displacement)

    if len(r_indexes) == 1:
#    line scan
        poscars = lineScan(poscar, displacement, nstep, r_indexes, m_indexes, g_indexes)
    elif len(r_indexes) == 2:
#    angle scan
        poscars = angleScan(poscar, displacement, nstep, r_indexes, m_indexes, g_indexes)
    elif len(r_indexes) == 3:
#    dihedral angle scan
        poscars = dihedralScan(poscar, displacement, nstep, r_indexes, m_indexes, g_indexes)
    else:
        print "input error of reference atom"
        sys.exit(2)

    makeScanJob(poscars)


if __name__ == "__main__":
    main()
#    file = "iPOSCAR"
#    poscar = POSCAR(file)

#    dihedral angle split
#   r: reference atom
#   m: motion atom
#   m_angle: motion angle
#   nstep: number of step
#   input
#    ref_indexes = [0, 1, 2]
#    mot_indexes = [3]
#    grp_indexes = [4, 5]
#    nstep = 8
#    angle = 90.0
#

