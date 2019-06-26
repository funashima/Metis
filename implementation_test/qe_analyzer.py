#!/usr/bin/env python3

from Metis.Espresso.QE2Spg import QE2Spg
import sys


if __name__ == '__main__':
    argc = len(sys.argv)
    if argc > 1:
        filename = sys.argv[1]
    else:
        print('usage QE2Spg.py [pw.x_file]')
        exit()
    spg = QE2Spg(filename).space_group
    spg.show_info()
