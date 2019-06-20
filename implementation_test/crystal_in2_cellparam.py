#!/usr/bin/env python3
#
# modified ParseConfig
#   written by Hiroki Funashima in Kobe 4 June 2019
#

from Metis.Base.ParseConfig import ParseConfig
import sys


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage {} [inputfile]'.format(sys.argv[0]))
        exit()
    ParseConfig(sys.argv[1])
