#!/usr/bin/env python3
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
import sys

def get_args():
    if len(sys.argv) < 2:
        print('usage {0} [ispg] [ichoice]'.format(sys.argv[0]))
        exit()
    ispg = sys.argv[1]
    if len(sys.argv) == 2:
        ichoice = 1
    else:
        ichoice = int(sys.argv[2])
    return [ispg, ichoice]


if __name__ == '__main__':
    ispg, ichoice = get_args()
    obj = GenerateSpaceGroup(ispg, ichoice)
    group_elements = obj.group_elements
    print()
    print('crystal system : {}'.format(obj.crystal_system))
    print('bravais lattice: {}'.format(obj.bravais_lattice_name))
    print('point group    : {}'.format(obj.point_group_name))
    print()
    print('Hermann-Maugin symbol:{}'.format(obj.hmname))
    print('Schoenfies symbol    :{}'.format(obj.schname))
    print()
    print('num of group elements = {}'.format(len(group_elements)))
    print('Group Elements')
    obj.display_group_elements()
    print()
    print()
    print('Group Table')
    print()
    obj.display_group_table()
