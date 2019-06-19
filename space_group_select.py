#!/home/funashima/local/bin/python3
from Metis.SpaceGroup.SearchSpaceGroup import SearchSpaceGroup
import os


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

min_spg_index = config_obj.min_spg_index
max_spg_index = config_obj.max_spg_index
natoms = config_obj.natoms

spg_obj = SearchSpaceGroup()

for ispg in range(min_spg_index, max_spg_index+1):
    print('Space Group #{}'.format(ispg))
    spg_obj.set_space_group(ispg)

    print(' -- possibility(atomic site) --')
    for wyckoff_patterns in spg_obj.get_pattern(natoms):
        print('  ', end='')
        print(wyckoff_patterns)

    print(' -- actual --')
    for wyckoff_patterns in spg_obj.set_atomic_position(natoms):
        for composition in wyckoff_patterns:
            print('  ', end='')
            for atomic_site in composition:
                wyckoff_name = spg_obj.get_wyckoff_name(atomic_site)
                print('{} '.format(wyckoff_name), end='')
            print()
    print()
