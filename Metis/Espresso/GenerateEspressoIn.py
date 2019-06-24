#!/usr/bin/env python3


class GenerateEspressoIn(object):
    def __init__(self, configfile,
                       ispg, ichoice=1,
                       atom_info):
        self.ispg = ispg
        self.ichoice = ichoice

        ParseConfigGenerateQEInput(configfile)

        self.get_crystal_structure():
        self.il = self.crystal_structure.il
        self.crystal_system = self.crystal_structure.crystal_system

    def get_crystal_structure(self):
        self.crystal_structure = CrystalStructure(self.ispg,
                                                  self.ichoice
                                                  max_coa_ratio=max_coa_ratio,
                                                  apf=apf, delta_apf=delta_apf,
                                                  thr_bond_length=thr_bond_length,
                                                  max_try=max_try,

    def set_calculation(self):
        pass

