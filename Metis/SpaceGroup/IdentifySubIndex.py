#!/usr/bin/env python3
from Metis.SpaceGroup.GenerateWyckoffPositionsList \
        import GenerateWyckoffPositionsList


class IdentifySubIndex(object):
    def __init__(self, natoms=None, ispg=None, target_list=None):
        self.natoms = natoms
        self.ispg = ispg
        self.target_list = target_list
        self.get_sub_index()

    def length_check(self, reagant, target_list):
        ex_reagant = len(self.decomp_wylist(reagant))
        ex_target_list = len(self.decomp_wylist(target_list))
        if ex_reagant % self.natoms != 0:
            print('====== Error(IdentifySubIndex) ======')
            print('regant:{} is inccorect.'.format(reagant))
            print('natoms = {}'.format(self.natoms))
            exit()
        if ex_target_list % self.natoms != 0:
            print('====== Error(IdentifySubIndex) ======')
            print('target_list:{} is inccorect.'.format(target_list))
            exit()
        if ex_reagant % ex_target_list != 0:
            print('====== Error(IdentifySubIndex) ======')
            print('incompatible.')
            print('target_list:{} is inccorect.'.format(target_list))
            exit()
        #
        # ex_reagant // ex_target_list: n of sub lattice
        #   e.g. face centered lattice = 4
        #        body centered lattice = 2
        #        base centered lattice = 2
        #
        return ex_reagant // ex_target_list

    def decomp_wylist(self, target_list, nlat=1):
        new_list = []
        for x in target_list:
            n = list(x)
            letter = n.pop()
            n = int(''.join(n))
            for i in range(n * nlat):
                new_list.append(letter)
        return sorted(new_list)

    def compare_wyckoff_lists(self, target_list, wyckoff_dic):
        reagant0 = wyckoff_dic[0]
        #
        # target is primitive or conventinal cell
        # (both case is ok )
        #
        nlat = self.length_check(reagant0, target_list)

        wylist1 = self.decomp_wylist(target_list, nlat=nlat)
        for (i, reagant) in enumerate(wyckoff_dic):
            wylist2 = self.decomp_wylist(reagant, nlat=nlat)
            if wylist1 == wylist2:
                return i
        return None

    def get_sub_index(self):
        wyckoff_obj = GenerateWyckoffPositionsList(min_spg_index=self.ispg,
                                                   natoms=self.natoms,
                                                   use_progress_bar=False)
        wyckoff_dic = wyckoff_obj.wyckoff_list[0]['wyckoff_positions']
        self.sub_index = self.compare_wyckoff_lists(self.target_list,
                                                    wyckoff_dic)
        return self
