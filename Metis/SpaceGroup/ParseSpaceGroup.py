#!/usr/bin/env python3
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff


class ParseSpaceGroup(object):
    def __init__(self, generator_file=None, wyckoff_data=None):
        self.generator_obj = ParseGenerator(generator_file)
        self.wyckoff_obj = ParseWyckoff(wyckoff_data)

    def set_space_group(self, ispg, ichoice=1):
        spg_obj = self.generator_obj.generator_list[ispg-1][ichoice-1]
        self.schname = spg_obj.schname
        self.hmname = spg_obj.hmname
        self.il = spg_obj.il
        wyckoff_obj = self.wyckoff_obj.wyckoff_position[ispg-1][ichoice-1]
        self.npos = wyckoff_obj.npos
        self.wyckoff_position_list = wyckoff_obj.position
        self.set_position_index_list()

    def check_freedom(self):
        self.wyckoff_freedom = {}
        for ipos in range(self.npos):
            atomic_site = self.wyckoff_position_list[ipos]['site_name']
            wyckoff_letter = atomic_site[-1]
            position = self.wyckoff_position_list[ipos]['position']
            self.wyckoff_freedom[wyckoff_letter] = False
            for x in position:
                if x in ['x', 'y', 'z', '-', 'w',
                         'd', 'm', 'p', 'n', 'r',
                         's', 'u', 'v']:
                    self.wyckoff_freedom[wyckoff_letter] = True

    def set_position_index_list(self):
        self.position_index_list = []
        self.position_dict = {}
        for ipos in range(self.npos):
            atomic_site = self.wyckoff_position_list[ipos]['site_name']
            wyckoff_letter = atomic_site[-1]
            n = int(atomic_site[:-1])
            if self.il == -1:  # rhombohedral lattice
                self.m = 3
            elif self.il == 0:  # hexagonal lattice
                self.m = 1
            elif self.il == 1:  # simple lattice
                self.m = 1
            elif self.il == 2:  # face centered lattice
                self.m = 4
            elif self.il == 3:  # bace centered lattice
                self.m = 2
            elif self.il == 4:  # c-center base centered lattice
                self.m = 2
            m = self.m
            if not n // m in self.position_index_list:
                self.position_index_list.append(n // m)
                self.position_dict[n//m] = []
            self.position_dict[n//m].append(wyckoff_letter)
            self.check_freedom()

    def get_letter2n(self, wyckoff_letter):
        for (n, letters) in self.position_dict.items():
            if wyckoff_letter in letters:
                return n * self.m
        return None

    def show_info(self):
        print('schname: {}'.format(self.schname))
        print('hmname: {}'.format(self.hmname))
        print('lattice type = {}'.format(self.il))
        for ipos in range(self.npos):
            wyckoff = self.wyckoff_position_list[ipos]
            print('  atmic_site:{}'.format(wyckoff['site_name']), end='')
            print('  {}'.format(wyckoff['position']))
        print('position index list = ', end='')
        print(self.position_index_list)
        print(self.position_dict)
        print(self.wyckoff_freedom)


if __name__ == '__main__':
    obj = ParseSpaceGroup(generator_file='generator',
                          wyckoff_data='wycoff')
    for ispg in range(1, 231):
        print('---')
        obj.set_space_group(ispg)
        obj.show_info()
