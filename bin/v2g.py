#!/usr/bin/env python
#
### modified date: 2013/11/21
#

import sys, os, getopt
from vtool.util import *

def main():
    def usage():
        print "Usage: v2g -i POSCAR [-o xxx.gjf] [-e H,C,O,...]"
        print " -h : help"
        print " -i : input file, ie POSCAR"
        print " -o : output file, ie xxx.gjf, xxx.com"
        print " -e : setup element type, ie H,C,O,..."

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], "hi:o:e:")

    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
 
    infile = None
    outfile = None
    elements = None
    for o, a in opt_list:
        if o in ('-h'):
            usage()
            sys.exit()
        elif o in ('-i'):
            infile = a
        elif o in ('-o'):
            outfile = a
        elif o in ('-e'):
            elements = a.split(',')
 
    if infile is None:
        print "No intput file"
        usage()
        sys.exit(2)

    p = POSCAR(infile)
    g = GJF()
    len(sys.argv)
    poscar2gjf(p, g, elements)

    if outfile is None:
        g.writeGJF()
    else:
        g.writeGJF(outfile)

if __name__ == "__main__":
    main()
