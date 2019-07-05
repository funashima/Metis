#!/usr/bin/env python3
from Metis.SpaceGroup.GenerateMultipleWyckoffPositionList \
        import GenerateMultipleWyckoffPositionList


if __name__ == '__main__':
    spg = 'Fd-3m'
    natom_list = [4, 8]
    wyckoff_list = GenerateMultipleWyckoffPositionList(ispg=spg,
                                                       natom_list=\
                                                       natom_list).\
        wyckoff_list
    for (i, wy) in enumerate(wyckoff_list):
        print('({0:>3d})  {1}'.format(i+1, wy))
