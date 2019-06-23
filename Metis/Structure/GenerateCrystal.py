#!/usr/bin/env python3
import math
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff
from Metis.Structure.GetRandomLatticeConstant import GetRandomLatticeConstant
from fractions import Fraction
import random
import pymatgen as mg


class GenerateCrystal(TspaceToolbox):
    def __init__(self, ispg=None, ichoice=1,
                 max_coa_ratio=2.0,
                 apf=0.52,
                 max_try=100,
                 atom_info=None):
        self.space_group = GenerateSpaceGroup(ispg=ispg, ichoice=ichoice)
        self.ispg = self.space_group.ispg
        self.ichoice = ichoice
        self.il = self.space_group.il
        self.max_coa_ratio = max_coa_ratio
        self.apf = apf
        self.max_try = max_try
        self.atom_info = atom_info
        self.get_atomic_positions()
        self.generate_full_atomic_positions()
        self.set_lattice()

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
        for atom in self.atom_info:
            atom_name = atom['element']
            natom = len(atom['positions'])
            r = self.get_radius(atom_name)
            volume += natom * self.calc_sphere_volume(r)
        if self.il in [0, 1]:
            volume *= 1
        elif self.il == 2:
            volume *= 4
        elif self.il in [3, 4]:
            volume *= 2
        elif self.il == -1:
            volume *= 3
        return volume / self.apf

    def set_lattice(self):
        const_volume = self.set_const_volume()
        crystal_system = self.space_group.crystal_system
        lattice = GetRandomLatticeConstant(const_volume=const_volume,
                                           crystal_system=crystal_system)
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
            wyckoff_letter_list = []
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
                atom_site = GenerateAtomicPosition(ispg=self.ispg,
                                                   ichoice=self.ichoice,
                                                   wyckoff_params=wyckoff_params)
                atomic_pos = atom_site.atomic_position
                atomic_pos_list.append(atomic_pos)
            self.atom_info[i]['positions'] = atomic_pos_list
        self.count_num_atoms()


    def count_num_atoms(self):
        for (j, atom) in enumerate(self.atom_info):
            element = atom['element']
            n = 0
            for i in range(len(atom['wyckoff_position'])):
                n += len(atom['positions'][i])
            self.atom_info[j]['natoms'] = n

    #
    #  test methods for basis
    #

    def test_basis(self):
        print('Basis set:')
        for i in range(3):
            print(' t[{}] = '.format(i+1), end='')
            for j in range(3):
                print(' {:9.6f}'.format(round(self.basis[i][j], 7)), end='')
            print()
        print()

        lat = ['a', 'b', 'c']
        ang = ['alpha', 'beta', 'gamma']
        print('--- self test for lattice constant ---')
        print(' reconstruct a, b, c, alpha, beta and gammma')
        for i in range(3):
            length = self.get_norm(self.basis[i])
            err = abs(length - self.lattice_length[i])
            print('  {0} = {1:8.6f} [ang] (err={2:5.3e})'.
                  format(lat[i], length, err))
        print('  c/a = {:8.5f}'.
              format(self.lattice_length[2] / self.lattice_length[0]))
        print()

        for i in range(3):
            j1 = (i+1) % 3
            j2 = (i+2) % 3
            angle = self.get_angle(self.basis[j1], self.basis[j2])
            err = abs(angle - self.lattice_angle[i])
            print('  {0:5s} = {1:8.5f} [deg] (err={2:5.3e})'.
                  format(ang[i], angle, err))
        print()
        print('  volume = {:10.6f} [ang^3]'.
              format(self.calc_unitcell_volume()))
        print()

    def test_atomic_position(self):
        print('-- Atomic Position(not including sublattice):')
        iatom = 0
        for atom in self.atom_info:
            element = atom['element']
            for (i, wyckoff_position) in enumerate(atom['wyckoff_position']):
                print('Atom:{}  wyckoff_position:{}'.
                      format(element, wyckoff_position))
                for atomic_position in atom['positions'][i]:
                    iatom += 1
                    print(' ({:>3d})'.format(iatom), end='')
                    for i in range(3):
                        print('  {:9.6f}'.
                              format(float(atomic_position[i])),
                              end='')
                    print()
        print()
        print('-- Atomic Position(including sublattice):')
        print('sublattice points in this lattice:')
        j = 0
        for sublattice in self.sub_lattice_pattern:
            j += 1
            print('  ({})'.format(j), end='')
            for i in range(3):
                print(' {:9.5f}'.format(float(sublattice[i])), end='')
            print()
        iatom = 0
        print('atomic coordinate:')
        for atom in self.full_atomic_position:
            iatom += 1
            element = atom['element']
            wyckoff_position = atom['wyckoff_position']
            atomic_position = atom['positions']
            print(' ({:>3d})'.format(iatom), end='')
            print(' atom:{:2s}'.format(element), end='')
            for i in range(3):
                print('  {:9.6f}'.
                      format(float(atomic_position[i])),
                      end='')
            print()
        print()

    def show_space_group_info(self):
        print('---- Space Group infomation ----')
        print('Crystal system :{}'.
              format(self.space_group.crystal_system))
        print('Bravais Lattice:{}'.
              format(self.space_group.bravais_lattice_name))
        print('Point Group:{}'.
              format(self.space_group.point_group_name))
        print()
        print('Index of Space group    :#{0} ({1})'.
              format(self.space_group.ispg,
                     self.space_group.ichoice))
        print('Schoenflies symbol      :{}'.
              format(self.space_group.schname))
        print('Hermann-Mauguin symbol  :{}'.
              format(self.space_group.hmname))
        print()
        print('Group Elements')
        self.space_group.display_group_elements()
        print()
        print('Group Table')
        self.space_group.display_group_table()
        print()

    def show_info(self):
        self.show_space_group_info()
        self.test_basis()
        self.test_atomic_position()
        self.show_neighbor_table()

    def show_neighbor_table(self):
        print('---- Neighbor Table ----')
        natom = len(self.full_atomic_position)
        for i in range(1, natom):
            atom_name0 = self.full_atomic_position[i-1]['element']
            print('From atom({0:>2d}):{1:>2s} '.format(i, atom_name0))
            for j in range(i+1, natom+1):
                atom_name1 = self.full_atomic_position[j-1]['element']
                length1 = self.calc_length(i, j, aunit=True)
                length0 = self.calc_length(i, j, aunit=False)
                print('  to atom({0:>2d}): {1:>2s}'.
                      format(j, atom_name1), end='')
                print(' length = {0:8.5f} ang'.format(length0), end='')
                print(' (= {0:6.5f})'.format(length1), end='')
                # sum of atomic radius
                r_atoms = 0
                for atom in [atom_name0, atom_name1]:
                    r_atoms = self.get_radius(atom)
                print('; sum of 2 atomic radii = {:7.5f} ang'.
                      format(r_atoms), end='')
                print('  state ... ', end='')
                diff = length0 - r_atoms
                ratio = (length0 / r_atoms) * 100
                if r_atoms >= length0:
                    print('NG', end='')
                else:
                    print('OK', end='')
                print(' ({0:8.5f} ang, {1:6.2f}%)'.format(diff, ratio))

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
