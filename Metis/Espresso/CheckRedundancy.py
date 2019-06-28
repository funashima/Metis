#!/usr/bin/env python3
import os
from Metis.Espresso.QE2Spg import QE2Spg


class CheckRedundancy(object):
    #
    # return value: [consistent_spg, space_group_name]
    #

    def check(self):
        self.prefix = '{0}_{1}_{2}'.format(self.compound_name,
                                           self.ispg,
                                           self.sub_index)
        if not os.path.isdir(self.prefix):
            print('directry:{} is not found.'.format(self.prefix))
            exit()
        self.analyze_true_symmetry()
        if self.full_symmetry['ispg'] == self.ispg:
            return [True, self.ispg]
        return [False, self.full_symmetry['ispg']]

    def analyze_true_symmetry(self):
        atom_info = []
        inputfile = self.prefix + '/espresso_relax.in'
        if not os.path.isfile(inputfile):
            print('===== Error(CheckRedundancy) =====')
            print('file:{} is not found.'.format(inputfile))
            exit()
        spg = QE2Spg(inputfile).space_group
        ispg = spg.ispg
        atoms = spg.atoms
        wyckoffs = spg.spg.get_symmetry_dataset()['wyckoffs']
        natoms = len(atoms)
        info = {}
        for i in range(natoms):
            atom_name = atoms[i]
            if atom_name not in info:
                info[atom_name] = []
            wyletter = wyckoffs[i]
            if wyletter not in info[atom_name]:
                info[atom_name].append(wyletter)
        atom_info = []
        for (element, wyckoff_position) in info.items():
            atom = {'element': element,
                    'wyckoff_position': wyckoff_position}
            atom_info.append(atom)

        self.full_symmetry = {'ispg': ispg, 'atom_info': atom_info}
        tempfile = 'espresso_relax.crystal_structure.in'
        if self.ispg == 221:
            exit()
        if os.path.isfile(tempfile):
            os.remove(tempfile)
        return self
