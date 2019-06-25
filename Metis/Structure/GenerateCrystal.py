#!/usr/bin/env python3
import math
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff
from Metis.Structure.GetRandomLatticeConstant import GetRandomLatticeConstant
from fractions import Fraction
import os
import random
import sys
import pymatgen as mg


class GenerateCrystal(TspaceToolbox):
    def __init__(self, ispg=None, ichoice=1,
                 max_coa_ratio=2.0,
                 apf=0.52,
                 delta_apf=0.9,
                 thr_bond_ratio=0.75,
                 max_try=100,
                 atom_info=None,
                 progress=True):
        self.space_group = GenerateSpaceGroup(ispg=ispg, ichoice=ichoice)
        self.ispg = self.space_group.ispg
        self.ichoice = ichoice
        self.il = self.space_group.il
        self.max_coa_ratio = max_coa_ratio
        self.apf = apf
        self.delta_apf = delta_apf
        self.max_try = max_try
        self.atom_info = atom_info
        self.thr_bond_ratio = thr_bond_ratio
        self.progress = progress
        self.main()

    def main(self):
        lattice_check = False
        itry = 0  # total count for trial
        max_jtry = 10  # max num of same apf
        jtry = 0  # reduce apf
        while not lattice_check:
            itry += 1
            if self.progress:
                sys.stdout.write('\r')
                sys.stdout.write('  tring to generate crystal # = {:>3d}'.format(itry))
                sys.stdout.write(' apf = {:5.3f}'.format(self.apf))
                sys.stdout.flush()
            jtry += 1
            if itry > self.max_try:
                break
            self.get_atomic_positions()
            self.generate_full_atomic_positions()
            self.set_lattice()
            lattice_check = self.check_reasonable_lattice()
            if jtry == max_jtry:
                jtry = 0
                self.apf *= self.delta_apf
        if self.progress:
            print()

    def calc_sphere_volume(self, r):
        return (4.0 * math.pi / 3.0) * (r**3)

    def get_radius(self, atomic_symbol):
        return float(mg.Element(atomic_symbol).atomic_radius)

    def set_const_volume(self):
        #
        # -Atomic Packing Factor(APF)-
        #
        # APF: bcc = pi*sqrt(3.0) / 8 = 0.68
        #      fcc/hcp = pi/sqrt(18) = 0.74
        #      diamond = 0.34
        #      simple cubic lattice = 0.52
        #
        volume = 0.0
        eps = 1.0e-8
        for atom in self.atom_info:
            atom_name = atom['element']
            natoms = atom['natoms']
            r = self.get_radius(atom_name)
            volume += natoms * self.calc_sphere_volume(r)
        if self.il in [0, 1]:
            volume *= 1
        elif self.il == 2:
            volume *= 4
        elif self.il in [3, 4]:
            volume *= 2
        elif self.il == -1:
            volume *= 3
        if abs(volume) < eps:
            print('Error: volume is too small !!')
            print('  volume = {}'.format(volume))
            exit()
        return volume / self.apf

    def set_lattice(self):
        const_volume = self.set_const_volume()
        self.crystal_system = self.space_group.crystal_system
        lattice = GetRandomLatticeConstant(const_volume=const_volume,
                                           crystal_system=self.crystal_system)
        a, b, c = lattice.get_lattice_length(max_coa_ratio=self.max_coa_ratio)
        alpha, beta, gamma = lattice.get_lattice_angle()
        self.lattice_length = [a, b, c]
        self.lattice_angle = [alpha, beta, gamma]
        self.generate_basis_vectors()

    def get_atomic_positions(self):
        wyckoff_letters = []
        atom_name_list = []
        for atom in self.atom_info:
            atom_name = atom['element']
            atom_name_list.append(atom_name)
            wyckoff_letters.append(atom['wyckoff_position'])

        nka = len(atom_name_list)
        wyckoff = ParseWyckoff()
        for i in range(nka):
            atomic_pos_list = []
            for wletter in wyckoff_letters[i]:
                wyckoff_params = []
                pos = wyckoff.get_atomic_position(ispg=self.ispg,
                                                  ichoice=self.ichoice,
                                                  wyckoff_letter=wletter)
                if any(['x' in x for x in pos]):
                    x = random.random()
                else:
                    x = None
                if any(['y' in x for x in pos]):
                    y = random.random()
                else:
                    y = None
                if any(['z' in x for x in pos]):
                    z = random.random()
                else:
                    z = None
                wyckoff_params = {'letter': wletter, 'variables': [x, y, z]}
                atom_site =\
                    GenerateAtomicPosition(ispg=self.ispg,
                                           ichoice=self.ichoice,
                                           wyckoff_params=wyckoff_params)
                atomic_pos = atom_site.atomic_position
                atomic_pos_list.append(atomic_pos)
            self.atom_info[i]['positions'] = atomic_pos_list
        self.count_num_atoms()

    def count_num_atoms(self):
        self.compound_name = ''
        for (j, atom) in enumerate(self.atom_info):
            atom_name = atom['element']
            n = 0
            for atom_pos in atom['positions']:
                n += len(atom_pos)
            self.atom_info[j]['natoms'] = n
            self.compound_name += '{0}{1}'.format(atom_name, n)

    #
    #  test methods for basis
    #

    def show_basis(self, filename=None):
        def sub_show_basis(fout):
            fout.write('Basis set:\n')
            for i in range(3):
                fout.write(' t[{}] = '.format(i+1))
                for j in range(3):
                    fout.write(' {:9.6f}'.format(round(self.basis[i][j], 7)))
                fout.write('\n')
            fout.write('\n')

            lat = ['a', 'b', 'c']
            ang = ['alpha', 'beta', 'gamma']
            fout.write('--- self test for lattice constant ---\n')
            fout.write(' reconstruct a, b, c, alpha, beta and gammma\n')
            for i in range(3):
                length = self.get_norm(self.basis[i])
                err = abs(length - self.lattice_length[i])
                fout.write('  {0} = {1:8.6f} [ang] (err={2:5.3e})\n'.
                           format(lat[i], length, err))
            fout.write('  c/a = {:8.5f}\n'.
                       format(self.lattice_length[2] / self.lattice_length[0]))
            fout.write('\n')

            for i in range(3):
                j1 = (i+1) % 3
                j2 = (i+2) % 3
                angle = self.get_angle(self.basis[j1], self.basis[j2])
                err = abs(angle - self.lattice_angle[i])
                fout.write('  {0:5s} = {1:8.5f} [deg] (err={2:5.3e})\n'.
                           format(ang[i], angle, err))
            fout.write('\n')
            fout.write('  volume = {:10.6f} [ang^3]\n'.
                       format(self.calc_unitcell_volume()))
            fout.write('\n')

        if filename is None:
            sub_show_basis(sys.stdout)
        else:
            with open(filename, 'a') as fout:
                sub_show_basis(fout)

    def show_atomic_position(self, filename=None):
        def sub_show_atomic_position(fout):
            if self.il > -1:
                fout.write('-- Atomic Position(not including sublattice):\n')
            else:
                fout.write('-- Atomic Position(not including sublattice):\n')
                fout.write('    (Hexagonal axis Coordinate)\n')
            iatom = 0
            for atom in self.atom_info:
                element = atom['element']
                for (i, wyckoff_position) in \
                        enumerate(atom['wyckoff_position']):
                    fout.write('Atom:{}  wyckoff_position:{}\n'.
                               format(element, wyckoff_position))
                    for atomic_position in atom['positions'][i]:
                        iatom += 1
                        fout.write(' ({:>3d})'.format(iatom))
                        for j in range(3):
                            fout.write('  {:9.6f}'.
                                       format(float(atomic_position[j])))
                        fout.write('\n')
            fout.write('\n')

            #
            # trigonal coordinate
            #   by Hiroki Funashima, 25 June 2019 in Kobe
            #
            a, _, c = self.lattice_length
            a_trg, alpha_trg = self.hex2trig_lattice_params(a, c)
            fout.write('\n')
            fout.write('Lattice constant in trigonal axis coordinate\n')
            fout.write('  a_trg = {:8.6f} [ang]\n'.format(a_trg))
            fout.write('  alpha_trg = {:8.6f} [deg]\n'.format(alpha_trg))
            jatom = 0
            if self.il == -1:
                fout.write('\n')
                fout.write('-- Atomic Position(primitive unit cell):\n')
                fout.write('    (trigonal axis Coordinate)\n')
                for atom in self.atom_info:
                    element = atom['element']
                    for (i, wyckoff_position) in \
                            enumerate(atom['wyckoff_position']):
                        fout.write('Atom:{}  wyckoff_position:{}\n'.
                                   format(element, wyckoff_position))
                        for atomic_position in atom['positions'][i]:
                            jatom += 1
                            fout.write(' ({:>3d})'.format(jatom))
                            trigonal = self.hex2trig(atomic_position)
                            for j in range(3):
                                fout.write('  {:9.6f}'.
                                           format(float(trigonal[j])))
                            fout.write('\n')
                fout.write('\n\n')
            fout.write('-- Atomic Position(including sublattice):\n')
            fout.write('sublattice points in this lattice:\n')
            j = 0
            for sublattice in self.sub_lattice_pattern:
                j += 1
                fout.write('  ({})'.format(j))
                for i in range(3):
                    fout.write(' {:9.5f}'.format(float(sublattice[i])))
                fout.write('\n')
            iatom = 0
            fout.write('atomic coordinate:\n')
            for atom in self.full_atomic_position:
                iatom += 1
                element = atom['element']
                wyckoff_position = atom['wyckoff_position']
                atomic_position = atom['positions']
                fout.write(' ({:>3d})'.format(iatom))
                fout.write(' atom:{:2s}'.format(element))
                for i in range(3):
                    fout.write('  {:9.6f}'.
                               format(float(atomic_position[i])))
                fout.write('\n')
            fout.write('\n')

        if filename is None:
            sub_show_atomic_position(sys.stdout)
        else:
            with open(filename, 'a') as fout:
                sub_show_atomic_position(fout)

    def show_space_group_info(self, filename=None):
        def sub_show_space_group_info(fout, filename):
            fout.write('---- Space Group infomation ----\n')
            fout.write('Crystal system :{}\n'.
                       format(self.space_group.crystal_system))
            fout.write('Bravais Lattice:{}\n'.
                       format(self.space_group.bravais_lattice_name))
            fout.write('Point Group:{}\n'.
                       format(self.space_group.point_group_name))
            fout.write('\n')
            fout.write('Index of Space group    :#{0} ({1})\n'.
                       format(self.space_group.ispg,
                              self.space_group.ichoice))
            fout.write('Schoenflies symbol      :{}\n'.
                       format(self.space_group.schname))
            fout.write('Hermann-Mauguin symbol  :{}\n'.
                       format(self.space_group.hmname))
            fout.write('\n')

        if filename is None:
            sub_show_space_group_info(sys.stdout, filename)
        else:
            with open(filename, 'w') as fout:
                sub_show_space_group_info(fout, filename)
        self.space_group.display_group_elements(filename=filename)
        self.space_group.display_group_table(filename=filename)

    def check_reasonable_lattice(self):
        natom = len(self.full_atomic_position)
        for i in range(1, natom):
            atom_name0 = self.full_atomic_position[i-1]['element']
            for j in range(i+1, natom+1):
                atom_name1 = self.full_atomic_position[j-1]['element']
                length = self.calc_length(i, j, aunit=False)
                r_atoms = 0
                for atom in [atom_name0, atom_name1]:
                    r_atoms = self.get_radius(atom)
                ratio = (length / r_atoms)
                if ratio < self.thr_bond_ratio:
                    return False
        return True

    def show_info(self, filename=None):
        if filename is not None:
            if os.path.isfile(filename):
                os.remove(filename)
        self.show_space_group_info(filename)
        self.show_basis(filename)
        self.show_atomic_position(filename)
        self.show_neighbor_table(filename)

    def show_neighbor_table(self, filename=None):
        def sub_show_neighbor_table(fout):
            fout.write('---- Neighbor Table ----\n')
            natom = len(self.full_atomic_position)
            for i in range(1, natom):
                atom_name0 = self.full_atomic_position[i-1]['element']
                fout.write('From atom({0:>2d}):{1:>2s}\n'.
                           format(i, atom_name0))
                for j in range(i+1, natom+1):
                    atom_name1 = self.full_atomic_position[j-1]['element']
                    length1 = self.calc_length(i, j, aunit=True)
                    length0 = self.calc_length(i, j, aunit=False)
                    fout.write('  to atom({0:>2d}): {1:>2s}'.
                               format(j, atom_name1))
                    fout.write(' length = {0:8.5f} ang'.format(length0))
                    fout.write(' (= {0:6.5f})'.format(length1))
                    # sum of atomic radius
                    r_atoms = 0
                    for atom in [atom_name0, atom_name1]:
                        r_atoms = self.get_radius(atom)
                    fout.write('; sum of 2 atomic radii = {:7.5f} ang'.
                               format(r_atoms))
                    fout.write('  state ... ')
                    diff = length0 - r_atoms
                    ratio = (length0 / r_atoms) * 100
                    if r_atoms >= length0:
                        fout.write('NG')
                    else:
                        fout.write('OK')
                    fout.write(' ({0:8.5f} ang, {1:6.2f}%)\n'.
                               format(diff, ratio))

        if filename is None:
            sub_show_neighbor_table(sys.stdout)
        else:
            with open(filename, 'a') as fout:
                sub_show_neighbor_table(fout)

    def calc_length(self, i, j, aunit=False):
        #
        # self.atomic_position[k] : wyckoff position
        # self.atom_list[k] : atom_type
        #
        u = self.full_atomic_position[i-1]['positions']
        v = self.full_atomic_position[j-1]['positions']
        if aunit:
            return self._measure_length(u, v) / self.lattice_length[0]
        else:
            return self._measure_length(u, v)

    def generate_full_atomic_positions(self):
        self.full_atomic_position = []
        self.sub_lattice_pattern = self.get_sub_lattices()

        for atom in self.atom_info:
            element = atom['element']
            for (i, wyckoff_position) in enumerate(atom['wyckoff_position']):
                for v in atom['positions'][i]:
                    info = {'element': element,
                            'wyckoff_position': wyckoff_position,
                            'positions': v}
                    self.full_atomic_position.append(info)
                    for sub_lattice in self.sub_lattice_pattern:
                        vv = self.add_sublattice(v, sub_lattice)
                        info = {'element': element,
                                'wyckoff_position': wyckoff_position,
                                'positions': vv}
                        self.full_atomic_position.append(info)


class GenerateAtomicPosition(TspaceToolbox):
    def __init__(self, ispg=None, ichoice=1, wyckoff_params=None):
        #
        # wyckoff params includes
        #   wyckoff letter and wyckoff params(x, y, z)
        #
        self.ichoice = ichoice
        self.space_group = GenerateSpaceGroup(ispg=ispg, ichoice=ichoice)
        self.ispg = self.space_group.ispg
        self.il = self.space_group.il
        self.atomic_position = []
        self.wyckoff_params = wyckoff_params
        self.wyckoff_obj = ParseWyckoff()
        self.generate_atomic_positions()

    def generate_atomic_positions(self):
        self.set_wyckoff_positions()
        for group_element in self.space_group.group_elements:
            for r in self.init_atoms:
                self.operate_space_group_op(r, group_element)

    def show_info(self):
        for atomic_position in self.atomic_position:
            print([round(float(x), 8) for x in atomic_position])

    def set_wyckoff_positions(self):
        self.init_atoms = []
        letter = self.wyckoff_params['letter']
        x, y, z = self.wyckoff_params['variables']
        pos = self.wyckoff_obj.get_atomic_position(ispg=self.ispg,
                                                   ichoice=self.ichoice,
                                                   x=x, y=y, z=z,
                                                   wyckoff_letter=letter)
        for i in range(3):
            if '/' in pos[i]:
                pos[i] = Fraction(pos[i])
            else:
                pos[i] = eval(pos[i])
        self.init_atoms.append(pos)

    def operate_space_group_op(self, r, group_element):
        rot = group_element['rot_matrix']
        tvec = group_element['translation']
        v = self.operate_spg_element(r, rot, tvec)
        if not self.check_existing_atomic_position(v):
            self.atomic_position.append(v)

    def check_existing_atomic_position(self, u):
        for position in self.atomic_position:
            if self.same_position(position, u):
                return True
        return False

    def same_position(self, u, v, eps=1.0e-8):
        if self.check_same_vectors(u, v, eps):
            return True
        for sub_lattice in self.get_sub_lattices():
            vv = self.add_sublattice(v, sub_lattice)
            if self.check_same_vectors(u, vv, eps):
                return True
        return False
