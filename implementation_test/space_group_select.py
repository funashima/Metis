#!/usr/bin/env python3
from Metis.SpaceGroup.SearchSpaceGroup import SearchSpaceGroup
import os
from tqdm import tqdm


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('file:{} is not found.'.format(configfile))
            exit()
        self.parse()

    def parse(self):
        self.min_spg_index = 1
        self.max_spg_index = 230
        for line in open(self.configfile, 'r'):
            linebuf = line.strip()
            if '#' in linebuf:
                linebuf.split('#')[0].strip()
            if linebuf == '':
                continue
            if '=' in linebuf:
                key, value = list(map(lambda x: x.strip(),
                                      linebuf.split('=')[:2]))
                key = key.lower()
                if key == 'min_spg_index':
                    self.min_spg_index = int(value)
                if key == 'max_spg_index':
                    self.max_spg_index = int(value)
                if key == 'natoms':
                    self.natoms = int(value)

config_obj = ParseConfig('space_group_select.in')


def get_wyckoff_list(space_group_list, natoms):
    wyckoff_list = []
    spg_obj = SearchSpaceGroup()
    for ispg in tqdm(space_group_list):
        spg_obj.set_space_group(ispg)
        pos = []
        for wyckoff_patterns in spg_obj.set_atomic_position(natoms):
            for composition in wyckoff_patterns:
                sub_pos = []
                for atomic_site in composition:
                    wyckoff_name = spg_obj.get_wyckoff_name(atomic_site)
                    sub_pos.append(wyckoff_name)
                pos.append(sub_pos)
        info = {'spg_index': ispg, 'wyckoff_positions': pos}
        wyckoff_list.append(info)
    return wyckoff_list


min_spg_index = config_obj.min_spg_index
max_spg_index = config_obj.max_spg_index
print('minimum index for space group = {}'.format(min_spg_index))
print('maximum index for space group = {}'.format(max_spg_index))

space_group_list = list(range(min_spg_index, max_spg_index + 1))

total_pattern = 0
print('now constructing database. please wait few seconds...')
for info in get_wyckoff_list(space_group_list, config_obj.natoms):
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
