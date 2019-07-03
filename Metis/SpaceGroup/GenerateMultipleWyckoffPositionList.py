#!/usr/bin/env python3
from Metis.SpaceGroup.ParseSpaceGroup import ParseSpaceGroup
from Metis.SpaceGroup.GenerateWyckoffPositionsList \
        import GenerateWyckoffPositionsList
import copy


class GenerateMultipleWyckoffPositionList(object):
    def __init__(self, ispg=None, natom_list=None):
        self.ispg = ispg
        self.natom_list = natom_list
        self.generator = ParseSpaceGroup().\
            set_space_group(ispg)
        self.main()

    def main(self):
        self.wyckoff_list = self.get_ex_wyckoff_list()

    def check_freedom(self, atomic_site):
        return self.generator.wyckoff_freedom[atomic_site[-1]]

    def get_wyckoff_list(self, natoms):
        return GenerateWyckoffPositionsList(min_spg_index=self.ispg,
                                            natoms=natoms).\
               wyckoff_list[0]['wyckoff_positions']

    def get_ex_wyckoff_list(self):
        if len(self.natom_list) == 1:
            return self.get_wyckoff_list(self.natom_list[0])
        wyckoff_list = []
        for i in range(1, len(self.natom_list)):
            if i == 1:
                gen_list = []
                for wy0 in self.get_wyckoff_list(self.natom_list[0]):
                    for wy1 in self.get_wyckoff_list(self.natom_list[1]):
                        total_list = [wy0, wy1]
                        flat_list = [x for row in total_list for x in row]
                        duplicate_list = [x for x in set(flat_list)
                                          if flat_list.count(x) > 1]
                        if duplicate_list:
                            has_not_freedom = \
                                    [not self.check_freedom(x) for
                                     x in duplicate_list]
                            if any(has_not_freedom):
                                continue
                        gen_list.append(total_list)
                if len(self.natom_list) == 2:
                    return gen_list
            else:
                if i > 2:
                    gen_list = copy.deepcopy(wyckoff_list)
                    wyckoff_list = []
                for wy0 in self.get_wyckoff_list(self.natom_list[i]):
                    for wy1 in gen_list:
                        total_wy1 = [x for row in wy1 for x in row]
                        flat_list = [x for row in [wy0, total_wy1]
                                     for x in row]
                        duplicate_list = [x for x in set(flat_list)
                                          if flat_list.count(x) > 1]
                        if duplicate_list:
                            has_not_freedom = [not self.check_freedom(x)
                                               for x in duplicate_list]
                            if any(has_not_freedom):
                                continue
                        new_comp = copy.deepcopy(wy1)
                        new_comp.append(wy0)
                        wyckoff_list.append(new_comp)
        return wyckoff_list
