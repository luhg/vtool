#!/usr/bin/env python
#
### modified date: 2013/12/19
#

import sys, os, getopt
from vtool.vasp import *

def main():
    def usage():
         print "Usage: vfreq.py -i OUTCAR [-o xxx.log]"
         print " -h : help"
         print " -i : input file, ie OUTCAR"
         print " -o : output file, ie xxx.log, xxx.out"

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

    o = OUTCAR(infile)
    o.writeLog(outfile)


if __name__ == "__main__":
    main()
