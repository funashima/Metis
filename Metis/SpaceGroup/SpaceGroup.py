#!/usr/bin/env python3
#
# ref: https://pymatgen.org/_modules/pymatgen/symmetry/analyzer.html
#      https://pymatgen.org/_modules/pymatgen/core/lattice.html
#

from Metis.Base.ParseConfig import ParseConfig
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class SpaceGroup(object):
    def __init__(self, configfile, symprec=0.01,  angle_trelance=6):
        self.config_obj = ParseConfig(configfile)
        self.symprec = symprec
        self.angle_trelance = angle_trelance
        self.main()

    def set_structure(self):
        #
        # ======= ATOMIC POSITIONS =======
        #
        a, b, c = self.config_obj.lattice_length[:3]
        alpha, beta, gamma = self.config_obj.lattice_angle[:3]
        lattice = mg.Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        self.atoms = self.config_obj.atom_list
        atomic_positions = self.config_obj.atomic_position
        self.structure = mg.Structure(lattice, self.atoms, atomic_positions)

    def main(self):
        self.set_structure()
        #
        # generate instance for space group
        #
        self.spg = SpacegroupAnalyzer(self.structure,
                                      symprec=self.symprec,
                                      angle_tolerance=self.angle_trelance)

    def show_info(self):
        print('----- symmetrized structure(conventional unit cell) -----')
        print(self.spg.get_symmetrized_structure())
        print()
        print('----- primitive structure -----')
        print(self.spg.find_primitive())
        print()
        #
        # name of space group
        #
        print('--- infomation of space group ---')
        print('  HM symbol: {} '.format(self.spg.get_space_group_symbol()))
        print('  Space Group number = #{}'.
              format(self.spg.get_space_group_number()))
        print('  point group : {}'.format(self.spg.get_point_group_symbol()))
        print()
        print()
        print('--- crystal system ---')
        print('  crystal system:{}'.format(self.spg.get_crystal_system()))
        print('  lattice type:{}'.format(self.spg.get_lattice_type()))

        #
        # total data set
        #
        dataset = self.spg.get_symmetry_dataset()
        wyckoff_position = dataset['wyckoffs']
        atom_pos = []
        atom_name = []
        wyckoff_data = []
        for (i, wyckoff) in enumerate(wyckoff_position):
            aname = self.atoms[i]
            if aname not in atom_name or wyckoff not in atom_pos:
                info = {'site_type_symbol': aname,
                        'wyckoff_letter': wyckoff,
                        'multiply': 0}
                wyckoff_data.append(info)
                atom_pos.append(wyckoff)
                atom_name.append(aname)
        for (i, wyckoff) in enumerate(wyckoff_position):
            aname = self.atoms[i]
            for info in wyckoff_data:
                if info['site_type_symbol'] == aname:
                    if info['wyckoff_letter'] == wyckoff:
                        info['multiply'] += 1
        print('  # of kinds of atoms = {}'.format(len(wyckoff_data)))
        print()
        for info in wyckoff_data:
            site_name = '{0}{1}-site'.format(info['multiply'],
                                             info['wyckoff_letter'])
            atom_name = info['site_type_symbol']
            print('  atom:{0}  wyckoff: {1}'.format(atom_name, site_name))
        return

    def symmetrized(self):
        self.get_symmetrized_structure()
