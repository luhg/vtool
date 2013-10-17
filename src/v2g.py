#!/usr/bin/env python
#
### modified date: 2013/07/26
#

from atm import *

if __name__ == "__main__":
    import sys
    import os

    file = sys.argv[1]
    p = POSCAR(file)
    g = GJF()

    len(sys.argv)
    if len(sys.argv) == 3:
        elements = sys.argv[2].split(',')
#        elements = ['A', 'B', 'C', 'D', 'E']
        poscar2gjf(p, g, elements)
    else:
        poscar2gjf(p, g)

    g.writeGJF()
