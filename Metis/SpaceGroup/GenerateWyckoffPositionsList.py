#!/usr/bin/env python3
#
# generation of wyckoff position list
#    written by Hiroki Funashima in Kobe, 2019
#
# this version is only 1 elements case...
#
from Metis.SpaceGroup.SearchSpaceGroup import SearchSpaceGroup
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
from tqdm import tqdm


class GenerateWyckoffPositionsList(object):
    def __init__(self, min_spg_index=8,
                 max_spg_index=None,
                 natoms=1,
                 use_progress_bar=False):
        spg = ParseGenerator()
        self.min_spg_index = spg.get_ispg(min_spg_index)

        if max_spg_index is None:
            self.max_spg_index = spg.get_ispg(min_spg_index)
        else:
            self.max_spg_index = spg.get_ispg(max_spg_index)
        self.natoms = natoms
        self.use_tqdm = use_progress_bar
        self.get_wyckoff_list()

    def get_wyckoff_list(self):
        self.wyckoff_list = []
        spg_obj = SearchSpaceGroup()
        space_group_list = list(range(self.min_spg_index,
                                      self.max_spg_index+1))
        if self.use_tqdm:
            spg_list = tqdm(space_group_list)
            print('now constructing database. please wait few seconds...')
            print()
        else:
            spg_list = space_group_list

        for ispg in spg_list:
            #
            # debug on 3 July 2019
            #
            #if 143 <= ispg <= 167:
            #    continue
            spg_obj.set_space_group(ispg)
            pos = []
            for wyckoff_patterns in spg_obj.set_atomic_position(self.natoms):
                for composition in wyckoff_patterns:
                    sub_pos = []
                    for atomic_site in composition:
                        wyckoff_name = spg_obj.get_wyckoff_name(atomic_site)
                        sub_pos.append(wyckoff_name)
                    pos.append(sub_pos)
            info = {'spg_index': ispg, 'wyckoff_positions': pos}
            self.wyckoff_list.append(info)
        return self

    def show_info(self):
        total_pattern = 0
        print()
        print('# of atoms = {}'.format(self.natoms))
        for info in self.wyckoff_list:
            ispg = info['spg_index']
            positions = info['wyckoff_positions']
            if len(positions) > 0:
                print('space group = #{}'.format(ispg))
                total_pattern += len(positions)
                for (i, pos) in enumerate(positions):
                    print(' ({0:>2d}) '.format(i+1), end='')
                    n = len(pos)
                    if n == 1:
                        print('{:>3s}'.format(pos[0]))
                    else:
                        for j in range(n):
                            print('{:>3s} '.format(pos[j]), end='')
                            if j < len(pos) - 1:
                                print('+ ', end='')
                            else:
                                print()
        print()
        print('-----')
        print('total {} patterns'.format(total_pattern))
        return self
