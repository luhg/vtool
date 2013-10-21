#!/usr/bin/env python
#
### modified date: 2013/07/26
#

import sys, os, getopt
from atm import *

def main():
    def usage():
        print "Usage: g2v -i xxx.gjf [-o POSCAR]"
        print " -h : help"
        print " -i : input file, ie xxx.gjf, xxx.com"
        print " -o : output file, ie POSCAR"

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], "hi:o:")

    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    infile = None
    outfile = None
    for o, a in opt_list:
        if o in ('-h'):
            usage()
            sys.exit()
        elif o in ('-i'):
            infile = a
        elif o in ('-o'):
            outfile = a

    if infile is None:
        print "No intput file"
        usage()
        sys.exit(2)

    g = GJF(infile)
    p = POSCAR()
    gjf2poscar(g, p)

#    print p.getLattice().getVectors()
    if outfile is None:
        p.writePOSCAR()
    else:
        p.writePOSCAR(outfile)

if __name__ == "__main__":
    main()
#    import sys
#    import os

#    file = sys.argv[1]
#    g = GJF(file)
#    p = POSCAR()
#    gjf2poscar(g, p)

#    print p.getLattice().getVectors()

#    p.writePOSCAR()

