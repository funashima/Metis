#!/usr/bin/env python3
#
# ref: https://pymatgen.org/_modules/pymatgen/symmetry/analyzer.html
#      https://pymatgen.org/_modules/pymatgen/core/lattice.html
#

from Metis.Base.ParseConfig import ParseConfig
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import re


class NotIdentifiedSpaceGroupError(Exception):
    ''' Exception Class, pymatgen cannot identify space group '''
    pass


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
        try:
            self.hmname = self.spg.get_space_group_symbol()
            self.ispg = self.spg.get_space_group_number()
        except TypeError:
            print()
            print(' **** INTERNAL LIBRARY WARNING ****')
            print('   pymatgen and spglib cannot identify space group', end='')
            print(' for smaller primitive cell.')
            print('   we expect this bug will be fixed near future.', end='')
            print('metis uses original primitive cell.')
            self.spg = None
            self.hmname = None
            self.ispg = None

    def generate_primitive_lattice(self):
        #
        # in this case, spglib cannot identify space group
        #
        if self.spg is None:
            return

        my_data = str(self.spg.find_primitive())
        data_region = False
        atomic_positions = []
        for line in my_data.split('\n'):
            linebuf = line.strip()
            if re.search('^abc\s+:', linebuf):
                lattice_length = \
                    [float(x) for x in linebuf.split(':')[1].split()]
            if re.search('^angles:', linebuf):
                lattice_angle = \
                    [float(x) for x in linebuf.split(':')[1].split()]
            if re.search('^#', linebuf) and 'SP' in linebuf:
                data_region = True
                continue
            if data_region:
                if re.search('^---', linebuf):
                    continue
                element = linebuf.split()[1]
                pos = [float(x) for x in linebuf.split()[2:]]
                info = {'element': element, 'position': pos}
                atomic_positions.append(info)
                lattice_type = self.hmname[0]
        compound_name = self.get_compound_name(atomic_positions)
        return {'ispg': self.ispg, 'hmname': self.hmname,
                'lattice_length': lattice_length,
                'lattice_angle': lattice_angle,
                'lattice_type': lattice_type,
                'atomic_positions': atomic_positions,
                'compound_name': compound_name}

    def get_compound_name(self, atomic_positions):
        natoms = None
        compound = ''
        pre_element = None
        for atom in atomic_positions:
            element = atom['element']
            if pre_element != element:
                if pre_element is not None:
                    compound += str(natoms)
                compound += element
                natoms = 1
                pre_element = element
            else:
                natoms += 1
        compound += str(natoms)
        return compound

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
        print('  HM symbol: {} '.format(self.hmname))
        print('  Space Group number = #{}'.format(self.ispg))
        print('  point group : {}'.format(self.spg.get_point_group_symbol()))
        print()
        print()
        print('--- crystal system ---')
        print('  crystal system:{}'.format(self.spg.get_crystal_system()))
        print('  lattice type:{}'.format(self.spg.get_lattice_type()))
        print()

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

    def get_conventional_cell(self):
        #
        # infomation of original primitive cell
        #
        self.prim_atom_info = {}
        data_set = self.spg.get_symmetry_dataset()
        wyckoffs = data_set['wyckoffs']
        natoms = len(self.atoms)
        atom_info_tmp = []
        for iatom in range(natoms):
            element = self.atoms[iatom],
            wyckoff_letter = wyckoffs[iatom]
            info = {'element': self.atoms[iatom],
                    'wyckoff_letter': wyckoffs[iatom]}
            atom_info_tmp.append(info)

        element_list = []
        atom_info = []
        for info in atom_info_tmp:
            element = info['element']
            wyckoff_letter = info['wyckoff_letter']
            if element not in element_list:
                element_list.append(element)
                info = {}
                info['element'] = element
                info['wyckoff_letter'] = [wyckoff_letter]
                atom_info.append(info)
            else:
                atom_info[-1]['wyckoff_letter'].append(wyckoff_letter)

        self.primitive_cell_info = {'ispg': self.ispg,
                                    'natoms': natoms,
                                    'atom_info': atom_info}

        #
        # conventional unit cell
        #
        data = str(self.spg.get_conventional_standard_structure())
        data_set = self.spg.get_symmetry_dataset()
        wyckoffs = data_set['wyckoffs']
        data_region = False
        atom_info = []
        self.atoms = []
        atomic_positions = []
        for line in data.split('\n'):
            linebuf = line.strip()
            if re.search('abc\s+:', linebuf):
                lattice_length = \
                    [float(x) for x in linebuf.split(':')[1].split()]
            if re.search('angles\s*:', linebuf):
                lattice_angle = \
                    [float(x) for x in linebuf.split(':')[1].split()]
            if re.search('^#\s+SP', linebuf):
                data_region = True
                continue
            if data_region:
                if re.search('^---', linebuf):
                    continue
                data_line = linebuf.split()
                element = data_line[1]
                position = [float(x) for x in data_line[2:]]
                self.atoms.append(element)
                atomic_positions.append(position)
        #
        # generate conventional infomation
        #
        a, b, c = lattice_length
        alpha, beta, gamma = lattice_angle
        lattice = mg.Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        try:
            self.structure = mg.Structure(lattice, self.atoms,
                                          atomic_positions)
            self.spg = SpacegroupAnalyzer(self.structure,
                                          symprec=self.symprec,
                                          angle_tolerance=self.angle_trelance)
            self.ispg = self.spg.get_space_group_number()
            self.hmname = self.spg.get_space_group_symbol()
        except TypeError:
            print()
            print(' **** INTERNAL LIBRARY WARNING ****')
            print('   pymatgen and spglib cannot identify space group', end='')
            print(' for smaller primitive cell.')
            print('   we expect this bug will be fixed near future.', end='')
            print(' metis uses original primitive cell.')
            self.spg = None
            self.ispg = None
            self.hmname = None
        return self

    def symmetrized(self):
        self.get_symmetrized_structure()
