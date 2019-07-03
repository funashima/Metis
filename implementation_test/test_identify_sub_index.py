#!/usr/bin/env python3
from Metis.SpaceGroup.GenerateWyckoffPositionsList \
        import GenerateWyckoffPositionsList
from Metis.SpaceGroup.IdentifySubIndex \
        import IdentifySubIndex
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
import random
import sys

min_spg = 8
max_spg = 230
max_natoms = 16
max_sub_index = 0
jspg = 0
hmname = ''
schname = ''

generator = ParseGenerator()

for natoms in range(1, max_natoms):
    total_pattern = 0
    for ispg in range(min_spg, max_spg+1):
        wyckoff_obj = GenerateWyckoffPositionsList(min_spg_index=ispg,
                                                   natoms=natoms,
                                                   use_progress_bar=False)
        wyckoff_dic = wyckoff_obj.wyckoff_list[0]['wyckoff_positions']
        for i in range(len(wyckoff_dic)):
            total_pattern += 1
            target_list = wyckoff_dic[i]
            random.shuffle(target_list)
            sys.stdout.write('\r')
            sys.stdout.write('> natoms = {0:>2d}  '.format(natoms))
            sys.stdout.write('ispg = {0:>3d} '.format(ispg))
            sys.stdout.write('sub_index = {0:>6d} '.format(i))
            if max_sub_index < i:
                max_sub_index = i
                jspg = ispg
                hmname = generator.get_hmname(ispg)
                schname = generator.get_schname(ispg)
            sys.stdout.write('(max of sub_index = {0:>6d} '.format(max_sub_index))
            sys.stdout.write('ispg = #{0:>3d}, '.format(jspg))
            sys.stdout.write('{:10s}'.format(hmname))
            sys.stdout.write('{:6s})'.format(schname))
            sys.stdout.write(' total {0:8d} patterns'.format(total_pattern))
            j = IdentifySubIndex(natoms=natoms,
                                 ispg=ispg,
                                 target_list=target_list).sub_index
            out_list = wyckoff_dic[j]
            if target_list != out_list:
                print('check is failed')
                print('natoms={}'.format(natoms))
                print('  ispg = {}'.format(ispg))
                print(i, target_list)
                print(j, out_list)
                exit()
    print()
print()
print('all test have been passed!')
