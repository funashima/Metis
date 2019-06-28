#!/usr/bin/env python3
from Metis.Structure.ParseConfigStructure import ParseConfigStructure
from Metis.SpaceGroup.SearchSpaceGroup import SearchSpaceGroup

class GetWyckoffList(object):
    def __init__(self, configfile):
        self.configure = ParseConfigStructure(configfile)
        print(self.configure.inputfile)

    def setup_wyckoff_list(self, min_spg_index=None,
                           max_spg_index=None):
        self.min_spg_index = min_spg_index
        self.max_spg_index = max_spg_index
        return self.expand_wyckoff_list()

    def get_wyckoff_list(self):
        self.wyckoff_list = []
        #
        # now only 1 atom cases...
        #
        natoms = self.configure.atom_info[0]['natoms']
        spg_obj = SearchSpaceGroup()
        min_spg_index = self.min_spg_index
        max_spg_index = self.max_spg_index
        for ispg in range(min_spg_index, max_spg_index+1):
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
            self.wyckoff_list.append(info)
        return self

    def expand_wyckoff_list(self):
        #  expand wyckoff_list to concurrent calc.
        self.get_wyckoff_list()
        search_pattern_lists = []
        for site_info in self.wyckoff_list:
            #
            # non only 1 atom cases...
            #
            element = self.configure.atom_info[0]['element']
            semi_core = self.configure.atom_info[0]['semi_core']
            ispg = site_info['spg_index']
            sub_index = 0
            for wyckoff_position in site_info['wyckoff_positions']:
                if(len(wyckoff_position) == 0):
                    break
                sub_index += 1
                atom_info = {'element': element,
                             'wyckoff_position': wyckoff_position,
                             'semi_core':semi_core}
                info = {'ispg': ispg, 'atom_info': atom_info,
                        'sub_index':sub_index,
                        'configfile':self.configfile}
                search_pattern_lists.append(info)
        return search_pattern_lists
