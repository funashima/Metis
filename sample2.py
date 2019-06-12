#!/usr/bin/env python3
from Metis.Structure.GetConstrainCondition import GetConstrainConditionCrystalIn
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage {0} inputfile'.format(sys.argv[0]))
        exit()
    filename = sys.argv[1]
    obj = GetConstrainConditionCrystalIn(filename,unit='bohr')
    print('volume = {}'.format(obj.volume))
    obj.check_bond_length()

