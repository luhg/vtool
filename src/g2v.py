#!/usr/bin/env python
#
### modified date: 2013/07/26
#

from atm import *

if __name__ == "__main__":
    import sys
    import os

    file = sys.argv[1]
    g = GJF(file)
    p = POSCAR()
    gjf2poscar(g, p)

#    print p.getLattice().getVectors()

    p.writePOSCAR()

