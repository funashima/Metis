#!/usr/bin/env python3
#
# sample demo file for SpacegroupAnalyzer
#   by Hiroki Funashima in Kobe, 2019
#
#
from Metis.SpaceGroup.SpaceGroup import SpaceGroup
import sys


if __name__ == '__main__':
    if len(sys.argv) > 1:
        configfile = sys.argv[1]
    else:
        configfile = 'crystal.in'
    SpaceGroup(configfile).show_info()
    #SpaceGroup(configfile).generate_cif()
