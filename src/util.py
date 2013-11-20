#!/usr/bin/env python
#
### modified date: 2013/11/21
#

from operator import itemgetter, attrgetter
from vtool.gaussian import *
from vtool.vasp import *

__author__ = ""
__date__ = "2013/11/21"

__version__ = "$Revision: 0.1$"


def gjf2poscar(gjf, poscar):
    """ Gaussain GJF convert VASP POSCAR
        gjf:    {GJF}
        poscar: {POSCAR}
    """
    poscar.setLattice(gjf._lattice_.getVectors() )
    for a in gjf._atoms_:
        poscar.addAtom(a)

#def poscar2gjf(poscar, gjf, elements = ['A', 'B', 'C', 'D']):
def poscar2gjf(poscar, gjf, elements = None):
    """ VASP POSCAR convert Gaussain GJF
        poscar:   {POSCAR}
        gjf:      {GJF}
        elements: {string array}
    """
    if elements != None:
        poscar.setElementsType(elements)
    for a in poscar._atoms_:
        gjf.addAtom(a)
    for v in poscar.getLattice().getVectors():
        tmp1, tmp2, tmp3 = v.getBasis()
        a = Atom('Tv', tmp1, tmp2, tmp3)
        gjf.addAtom(a)


if __name__ == "__main__":
    import sys
    import os

    pass
