#!/usr/bin/env python3
from Metis.Base.ParseConfig import ParseConfig
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage {} crystal_structure file'.format(sys.argv[0]))
        exit()

    filename = sys.argv[1]
    config_obj = ParseConfig(filename)
    print()
    natom = config_obj.natom
    print('num of atoms in unit cell = {}'.format(natom))
    print()
    for i in range(3):
        print('Basis:{}'.format(config_obj.basis[i]))
    for i in range(natom):
        print('atom:{}  {}'.format(config_obj.atom_list[i], config_obj.atomic_position[i]))
    print()
    for i in range(natom):
        print('{:2s} '.format(config_obj.atom_list[i]), end='')
        for j in range(i):
            print('                  ', end='')
        for j in range(i, natom):
            #print('(i,j) = ({0}, {1}), length = {2:5.3f}'.format(i, j, config_obj.calc_length(i, j)), end='')
            print('length = {0:8.5f} '.format(config_obj.calc_length(i, j, aunit=True)), end='')
        print()
