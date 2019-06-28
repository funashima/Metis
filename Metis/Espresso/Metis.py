#!/usr/bin/env python3
from Metis.Espresso.GenerateEspressoIn import GenerateEspressoIn
from Metis.SpaceGroup.GetWyckoffList import GetWyckoffList


class Metis(GetWyckoffList):
    def __init__(self, configfile):
        super().__init__(configfile)
        self.configfile = configfile
        self.nnode = self.configure.nnode
        self.submit_job = self.configure.submit_job
        self.main()

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
        min_spg_index = self.configure.min_spg_index
        max_spg_index = self.configure.max_spg_index
        search_pattern_lists =\
            self.setup_wyckoff_list(min_spg_index=min_spg_index,
                                    max_spg_index=max_spg_index,)

        for target in search_pattern_lists:
            self.generate_crystal(target)
