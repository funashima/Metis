#!/usr/bin/env python3
#
# cif2space group
#   written by Hiroki Funashima in Kobe, 29 May 2019
#

from Metis.SpaceGroup.Cif2Spg import Cif2Spg
from Metis.SpaceGroup.SpaceGroup import SpaceGroup
import os
import sys


if __name__ == '__main__':
    argc = len(sys.argv)
    if argc > 1:
        ciffile = sys.argv[1]
    else:
        print('usage Cif2Spg.py [cif_file]')
        exit()
    name, ext = os.path.splitext(os.path.basename(ciffile))
    configfile = name + '.crystal_structure.in'
    Cif2Spg(ciffile, outfile=configfile)
    SpaceGroup(configfile).show_info()
