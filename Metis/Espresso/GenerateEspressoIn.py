#!/usr/bin/env python3


class GenerateEspressoIn(object):
    def __init__(self, configfile,
                       ispg, ichoice=1,
                       atom_info):
        self.crystal_structure = CrystalStructure(ispg,
                                                  ichoice
                                                  max_coa_ratio=max_coa_ratio,
                                                  apf=apf, delta_apf=delta_apf,
                                                  thr_bond_length=thr_bond_length,
                                                  max_try=max_try,
                                                  atom_info=atom_info)
        self.il = self.crystal_structure.il
        self.crystal_system = self.crystal_structure.crystal_system

    def hogehoge(self):
        pass
