#!/usr/bin/env python3
#
# interface for GenerateWyckoffPositionList
#   for Metis written by Hiroki Funashima, 2019 in Kobe
#

from Metis.Structure.ParseConfigStructure import ParseConfigStructure
from Metis.SpaceGroup.GenerateWyckoffPositionsList \
        import GenerateWyckoffPositionsList


class GetWyckoffList(object):
    def __init__(self, configfile):
        self.configure = ParseConfigStructure(configfile)

    def setup_wyckoff_list(self, min_spg_index=None,
                           max_spg_index=None):
        self.min_spg_index = min_spg_index
        self.max_spg_index = max_spg_index
        return self.expand_wyckoff_list()

    def get_wyckoff_list(self):
        #
        # now only 1 atom cases...
        #
        natoms = self.configure.atom_info[0]['natoms']
        self.wyckoff_list = GenerateWyckoffPositionsList(
                            min_spg_index=self.min_spg_index,
                            max_spg_index=self.max_spg_index,
                            natoms=natoms,
                            use_progress_bar=False).wyckoff_list
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
                             'semi_core': semi_core}
                info = {'ispg': ispg, 'atom_info': atom_info,
                        'sub_index': sub_index,
                        'configfile': self.configfile}
                search_pattern_lists.append(info)
        return search_pattern_lists
