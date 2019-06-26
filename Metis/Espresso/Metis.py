#!/usr/bin/env python3
from concurrent import futures
from Metis.Structure.ParseConfigStructure import ParseConfigStructure
from Metis.SpaceGroup.SearchSpaceGroup import SearchSpaceGroup
from Metis.Espresso.GenerateEspressoIn import GenerateEspressoIn

class Metis(object):
    def __init__(self, configfile):
        self.configfile = configfile
        self.configure = ParseConfigStructure(configfile)
        self.nnode = self.configure.nnode
        self.submit_job = self.configure.submit_job
        self.main()

    def get_wyckoff_list(self):
        self.wyckoff_list = []
        #
        # non only 1 atom cases...
        #
        natoms = self.configure.atom_info[0]['natoms']
        spg_obj = SearchSpaceGroup()
        min_spg_index = self.configure.min_spg_index
        max_spg_index = self.configure.max_spg_index
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
        self.search_pattern_lists = []
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
                self.search_pattern_lists.append(info)

    def generate_crystal(self, pattern):
        ispg = pattern['ispg']
        atom_info = [pattern['atom_info']]
        sub_index = pattern['sub_index']
        configfile = pattern['configfile']
        GenerateEspressoIn(configfile=configfile,
                           qe_inputfile='espresso_relax.in',
                           ispg=ispg,
                           sub_index=sub_index,
                           atom_info=atom_info,
                           submit_job=self.submit_job)


    def main(self):
        ncpu = self.nnode
        self.expand_wyckoff_list()
        future_list = []
        for target in self.search_pattern_lists:
            self.generate_crystal(target)
        return

        with futures.ProcessPoolExecutor(max_workers=ncpu) as executer:
            for target in [self.search_pattern_lists for j in range(ncpu)]:
                future = executer.submit(fn=self.generate_crystal, pattern=target)
                future_list.append(future)
            _ = futures.as_completed(fs=future_list)
