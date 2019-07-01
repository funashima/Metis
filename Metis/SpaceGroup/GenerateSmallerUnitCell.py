#!/usr/bin/env python3
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.Espresso.QE2Spg import QE2Spg
from Metis.SpaceGroup.SpaceGroup import SpaceGroup
from Metis.SpaceGroup.GenerateWyckoffPositionsList import GenerateWyckoffPositionsList
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
        self.spg = SpaceGroup(filename).get_conventional_cell()
        if os.path.isfile(filename):
            os.remove(filename)
        self.spg.show_info()
        #
        #  identified reduced cell 
        #

        #
        # experimental!!
        #
        wyckoff = GenerateWyckoffPositionsList(min_spg_index = self.spg,
                                               max_spg_index = self.spg,
                                               natoms= self.spg.natoms)

