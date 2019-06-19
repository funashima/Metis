#!/bin/env python3
#
# nop: number of symmetry operations
#
from Metis.SpaceGroup.ParseSpaceGroup import ParseSpaceGroup


class SearchSpaceGroup(object):
    def __init__(self):
        self.spg_info = ParseSpaceGroup()

    def set_space_group(self, ispg, ichoice=1):
        self.spg_info.set_space_group(ispg)
        self.nop = self.spg_info.position_index_list
        self.wyckoff_freedom = self.spg_info.wyckoff_freedom
        self.show_spginfo()

    def get_pattern(self, n):
        composition = []
        for x in self._make_combination(n):
            if x not in composition:
                composition.append(x)
        return composition

    def show_spginfo(self):
        self.schname = self.spg_info.schname
        self.hmname = self.spg_info.hmname
        il = self.spg_info.il
        print('Shoenflies Symbol:{}'.format(self.schname))
        print('Hermann-Mauguin Symbol:{}'.format(self.hmname))
        if il == -1:
            lat = 'rhombohedlral'
        elif il == 0:
            lat = 'hexagonal'
        elif il == 1:
            lat = 'simple'
        elif il == 2:
            lat = 'face centered'
        elif il == 3:
            lat = 'body centered'
        elif il == 4:
            lat = 'c-center base centered'
        print('Lattice Type: {} lattice'.format(lat))

    def check_duplicate(self, xx):
        new_list = []
        for pattern in xx:
            uniq_list = list(set(pattern))
            black_list = []
            for letter in uniq_list:
                if not self.wyckoff_freedom[letter]:
                    black_list.append(letter)
            check = True
            for letter in black_list:
                if pattern.count(letter) > 1:
                    check = False
            if check:
                new_list.append(pattern)
        return new_list

    def set_atomic_position(self, n):
        composition = []
        for pattern in self.get_pattern(n):
            dict_pattern = []
            for atomic_site in pattern:
                x = self.spg_info.position_dict[atomic_site]
                dict_pattern.append(x)
            nn = len(dict_pattern)
            if nn == 1:
                xx = dict_pattern[0]
            else:
                xx = []
                for arya in dict_pattern[0]:
                    for aryb in dict_pattern[1]:
                        c = list(arya) + list(aryb)
                        c.sort()
                        if c not in xx:
                            xx.append(c)
                xx = self.check_duplicate(xx)
                if nn >= 3:
                    for j in range(2, nn):
                        yy = []
                        for arya in xx:
                            for aryb in dict_pattern[j]:
                                c = arya + list(aryb)
                                c.sort()
                                if c not in xx:
                                    yy.append(c)
                        xx = []
                        for y in self.check_duplicate(yy):
                            if y not in xx:
                                xx.append(y)
                xx = self.check_duplicate(xx)
            if len(xx) > 0:
                composition.append(xx)
        return composition

    def get_wyckoff_name(self, letter):
        n = self.spg_info.get_letter2n(letter)
        return str(n) + letter

    def _make_combination(self, n):
        combination_array = []
        if n == min(self.nop):
            combination_array = [[min(self.nop)]]
        else:
            for m in self.nop:
                if m < min(self.nop):
                    continue
                if n < m:
                    continue
                if n == m:
                    combination_array.append([m])
                else:
                    for ary in self._make_combination(n-m):
                        new_ary = [m] + ary
                        new_ary.sort()
                        combination_array.append(list(reversed(new_ary)))
        return combination_array
