#!/usr/bin/env python3
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.Espresso.QE2Spg import QE2Spg
import os
from fractions import Fraction


class GenerateSmallerUnitCell(TspaceToolbox):
    def __init__(self, filename):
        spg = QE2Spg(filename).space_group
        self.atom_info = spg.generate_primitive_lattice()
        self.remove_tmpfile()
        self.set_atom_data()

    def remove_tmpfile(self):
        tmpfile = 'espresso_relax.crystal_structure.in'
        if os.path.isfile(tmpfile):
            os.remove(tmpfile)

    def set_atom_data(self):
        #
        # coordinate is not conventional unit but primitive unit
        #
        lattice_length = self.atom_info['lattice_length']
        lattice_angle = self.atom_info['lattice_angle']
        self.lattice_type = self.atom_info['lattice_type']
        atomic_positions = self.atom_info['atomic_positions']
        # set ntyp and nat
        #self.check_natoms(atomic_positions)
        #self.primitive_to_conventional(atomic_positions,
        #                               lattice_length,
        #                               lattice_angle)
        self.generate_prim_crystal_in(lattice_length,
                                      lattice_angle,
                                      atomic_positions)

    def generate_prim_crystal_in(self,
                                 lattice_length,
                                 lattice_angle,
                                 atomic_positions,
                                 filename='prim_crystal.in'):
        with open(filename, 'w') as fout:
            fout.write('#\n')
            fout.write('# lattice constant\n')
            fout.write('#\n')
            fout.write('\n')
            fout.write('lattice_length = ')
            for (i, x) in enumerate(lattice_length):
                fout.write('{:10.5f}'.format(x))
                if i < 2:
                    fout.write(',')
                else:
                    fout.write('\n')
            fout.write('lattice_angle  = ')
            for (i, x) in enumerate(lattice_angle):
                fout.write('{:10.5f}'.format(x))
                if i < 2:
                    fout.write(',')
                else:
                    fout.write('\n')
            fout.write('\n')
            fout.write('#\n')
            fout.write('# atomic position in conventional unit cell\n')
            fout.write('#\n')
            fout.write('begin_atomlist:\n')
            for atom in atomic_positions:
                element = atom['element']
                fout.write('    atom_type = {:2s} ; '.format(element))
                fout.write('position =')
                for (i, x) in enumerate(atom['position']):
                    fout.write('{:11.7f}'.format(x))
                    if i < 2:
                        fout.write(',')
                    else:
                        fout.write('\n')
            fout.write('end_atomlist:\n')

    def check_natoms(self, atomic_positions):
        self.natoms = len(atomic_positions)
        atom_array = []
        for atom in atomic_positions:
            element = atom['element']
            if atom_array not in atom_array:
                atom_array.push(element)
        self.ntype = len(atom_array)

    def primtive_to_conventional(self, atomic_positions,
                                 lattice_length,
                                 lattice_angle):
        self.atomic_positions = []
        if self.lattice_type in ['P', 'R']:
            rot_matrix, inv_rot_matrix = self.to_conv_P()
        elif self.lattice_type == 'F':
            rot_matrix, inv_rot_matrix = self.to_conv_F()
        elif self.lattice_type == 'I':
            rot_matrix, inv_rot_matrix = self.to_conv_I()
        elif self.lattice_type == 'C':
            rot_matrix, inv_rot_matrix = self.to_conv_C()
        else:
            print('unknown lattice_type:{}'.format(self.lattice_type))
            exit()
        self.to_convert_base(atomic_positions, rot_matrix)
        self.calc_lattice_params(inv_rot_matrix, lattice_length, lattice_angle)

    def get_oblique_inner_product(self, inv_rot_matrix,
                                  lattice_length, lattice_angle):

        #
        # length
        #
        for i in range(3):
            self.lattice_length[i] = 0
            n = [inv_rot_matrix[0][i],
                 inv_rot_matrix[1][i],
                 inv_rot_matrix[2][i]]
            for j in range(3):
                if i == j:
                    cos_theta = 1
                else:
                    cos_theta = self.deg2cosine(lattice_angle[3-i-j])
                self.lattice_length[i] += (n[i] * lattice_length[i]) *\
                                          (n[j] * lattice_length[j]) *\
                    cos_theta

        #
        # angle
        #
        for i in range(3):
            pass

    def to_convert_base(self, atomic_positions, rot_matrix):
        self.atomic_positions = []
        for atom in atomic_positions:
            element = atom['element']
            r = atom['position']
            pos = self.matrix_dot_vector(rot_matrix, r)
            info = {'element': element, 'position': pos}
        self.atomic_positions.append(info)

    def to_convert_P(self):
        rot_matrix = [[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]]
        return [rot_matrix, rot_matrix]

    def to_convert_F(self):
        half = Fraction('1/2')
        rot_matrix = [[0, half, half],
                      [half, 0, half],
                      [half, half, 0]]
        inv_rot_matrix = [[-1, 1, 1],
                          [1, -1, 1],
                          [1, 1, -1]]
        return [rot_matrix, inv_rot_matrix]

    def to_convert_I(self):
        half = Fraction('1/2')
        rot_matrix = [[-half, half, half],
                      [half, -half, half],
                      [half, half, -half]]
        inv_rot_matrix = [[0, 1, 1],
                          [1, 0, 1],
                          [1, 1, 0]]
        return [rot_matrix, inv_rot_matrix]

    def to_convert_C(self):
        half = Fraction('1/2')
        rot_matrix = [[half, half, 0],
                      [half, -half, 0],
                      [0, 0, 1]]
        inv_rot_matrix = [[1, 1, 0],
                          [1, -1, 0],
                          [0, 0, 1]]
        return [rot_matrix, inv_rot_matrix]
