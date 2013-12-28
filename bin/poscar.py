#!/usr/bin/env python
#
### modified date: 2013/12/29
#

import sys, os, getopt
from vtool.util import *

def main():
    def usage():
        print "Usage: poscar -i CONTCAR [-o POSCAR]"
        print " -h : help"
        print " -i : input file, ie CONTCAR, POSCAR (Automatically detect Cartesian or Direct)"
        print " -o : output file, ie POSCAR"
        print " These should be overridden if your input and output file name are the same."

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

    p = POSCAR(infile)
    if p.getCoorndinateType()[0].upper() == 'D':
        p.directToCartesian()
    elif p.getCoorndinateType()[0].upper() == 'C':
        p.cartesianToDirect()
    else:
        pass

    if outfile is None:
    	p.writePOSCAR()
    else:
        p.writePOSCAR(outfile)

if __name__ == "__main__":
    main()
