#!/usr/bin/env python3
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.Espresso.QE2Spg import QE2Spg
from Metis.SpaceGroup.SpaceGroup import SpaceGroup
#from Metis.SpaceGroup.GenerateWyckoffPositionsList \
#        import GenerateWyckoffPositionsList
from Metis.SpaceGroup.IdentifySubIndex \
        import IdentifySubIndex
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
        #
        # generate crystal data (temporary file)
        #
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
        #
        # transform primitive cell to conventional cell
        #
        self.spg = SpaceGroup(filename).get_conventional_cell()
        if os.path.isfile(filename):
            os.remove(filename)
        
        #
        #  identified reduced cell 
        #

        #
        # experimental!!
        #
        natoms = self.spg.primitive_cell_info['natoms']
        # only 1 element ....
        self.atom_info = []
        for atom in self.spg.primitive_cell_info['atom_info']:
            info = {'element': atom['element'],
                    'natoms': len(atom['wyckoff_letter']),
                    'wyckoff_position': sorted(set(atom['wyckoff_letter']))}
            self.atom_info.append(info)

        for iatom in range(1):
            atom_info = self.spg.primitive_cell_info['atom_info'][iatom]
            target_list = atom_info['wyckoff_letter']
            element = atom_info['element']
            sub_index = IdentifySubIndex(natoms=natoms,
                                        ispg=self.spg.ispg,
                                        target_list=target_list).sub_index
            self.ispg = self.spg.ispg
            self.sub_index = sub_index+1
            self.dirname = '{0}{1}_{2}_{3}'.format(element, natoms, self.spg.ispg, sub_index+1)
            self.compound_name = '{0}{1}'.format(element, natoms)
