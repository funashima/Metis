#!/usr/bin/env python3
import math
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff
from Metis.Structure.GetRandomLatticeConstant import GetRandomLatticeConstant
from fractions import Fraction
import random
import pymatgen as mg


class GenerateCrystal(object):
    def __init__(self, ispg=None, ichoice=1,
                 max_coa_ratio=4.0,
                 atom_info=None):
        self.space_group = GenerateSpaceGroup(ispg=ispg, ichoice=ichoice)
        self.ispg = self.space_group.ispg
        self.ichoice = ichoice
        self.il = self.space_group.il
        self.max_coa_ratio = max_coa_ratio
        self.atom_info = atom_info
        self.get_atomic_positions()
        self.set_lattice()

    def calc_sphere_volume(self, r):
        return (4.0 * math.pi / 3.0) * (r**3)

    def set_const_volume(self):
        volume = 0.0
        for atom in self.atom_data:
            atom_name = atom['name']
            natom = len(atom['positions'])
            element = mg.Element(atom_name)
            r = element.atomic_radius
            volume += natom * self.calc_sphere_volume(r)
        if self.il in [0, 1]:
            volume *= 1
        elif self.il == 2:
            volume *= 4
        elif self.il in [3, 4]:
            volume *= 2
        elif self.il == -1:
            volume *= 3
        return volume

    def set_lattice(self):
        const_volume = self.set_const_volume()
        crystal_system = self.space_group.crystal_system
        lattice = GetRandomLatticeConstant(const_volume=const_volume,
                                           crystal_system=crystal_system)
        a, b, c = lattice.get_lattice_length(max_coa_ratio=self.max_coa_ratio)
        alpha, beta, gamma = lattice.get_lattice_angle()
        self.lattice_param = {'a': a, 'b': b, 'c': c,
                              'alpha': alpha, 'beta': beta, 'gamma': gamma}

    def get_atomic_positions(self):
        wyckoff_letters = []
        atom_name_list = []
        self.atom_data = []
        for atom in self.atom_info:
            atom_name = atom['name']
            atom_site = atom['wyckoff_position']
            atom_name_list.append(atom_name)
            wyckoff_letters.append(atom_site)

        nka = len(atom_name_list)
        wyckoff = ParseWyckoff()
        for i in range(nka):
            wletter = wyckoff_letters[i]
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
            self.atom_data.append({'name': atom_name[i],
                                   'wyckoff_position': atom_site,
                                   'positions': atomic_pos})


class GenerateAtomicPosition(object):
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
        print(self.init_atoms)
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
        half = Fraction('1/2')
        tr1 = Fraction('1/3')
        tr2 = 2 * tr1
        if self.check_same_vectors(u, v, eps):
            return True
        if self.il == 2:
            pattern = [[0, half, half], [half, 0, half], [half, half, 0]]
            for sub_lattice in pattern:
                vv = self.add_sublattice(v, sub_lattice)
                if self.check_same_vectors(u, vv, eps):
                    return True
        elif self.il == 3:
            sub_lattice = [half, half, half]
            vv = self.add_sublattice(v, sub_lattice)
            if self.check_same_vectors(u, vv, eps):
                return True
        elif self.il == 4:
            sub_lattice = [half, half, 0]
            vv = self.add_sublattice(v, sub_lattice)
            if self.check_same_vectors(u, vv, eps):
                return True
        elif self.il == -1:
            pattern = [[tr1, tr2, tr2], [tr2, tr1, tr1]]
            for sub_lattice in pattern:
                vv = self.add_sublattice(v, sub_lattice)
                if self.check_same_vectors(u, vv, eps):
                    return True

    def add_sublattice(self, u, sub_lattice):
        v = [0, 0, 0]
        for i in range(3):
            v[i] = u[i] + sub_lattice[i]
        return self.reduce_coordinate(v)

    def check_same_vectors(self, u, v, eps=1.0e-6):
        uu = self.reduce_coordinate(u)
        vv = self.reduce_coordinate(v)
        for i in range(3):
            if abs(uu[i] - vv[i]) > eps:
                return False
        return True

    def reduce_coordinate(self, u, eps=1.0e-6):
        half = Fraction('1/2')
        for i in range(3):
            while abs(u[i]) > half:
                if u[i] < -half:
                    u[i] += 1
                if u[i] > half:
                    u[i] -= 1
            if abs(u[i] + half) < eps:
                u[i] = half
        return u

    def tsp2frac(self, tvec):
        x = Fraction('{0}/{1}'.format(tvec[0], tvec[1]))
        y = Fraction('{0}/{1}'.format(tvec[2], tvec[3]))
        z = Fraction('{0}/{1}'.format(tvec[4], tvec[5]))
        return [x, y, z]

    def frac2tsp(self, frac):
        ary = [str(x).split('/') for x in frac]
        for (i, x) in enumerate(ary):
            if len(x) == 1:
                ary[i] = [0, 1]
            else:
                ary[i] = [int(y) for y in x]
        return [flatten for inner in ary for flatten in inner]

    def to_matrix_representation(self, code):
        if code.strip().lower() == 'x':
            return [1, 0, 0]
        if code.strip().lower() == 'y':
            return [0, 1, 0]
        if code.strip().lower() == 'z':
            return [0, 0, 1]
        if code.strip().lower() == 'w':
            return [1, -1, 0]
        if code.strip().lower() == '-x':
            return [-1, 0, 0]
        if code.strip().lower() == '-y':
            return [0, -1, 0]
        if code.strip().lower() == '-z':
            return [0, 0, -1]
        if code.strip().lower() == '-w':
            return [-1, 1, 0]
        self.unknown_operator

    def operate_spg_element(self, u, rot, t):
        #
        # -- in --
        # u: atomic position
        # rot: rotational matrix
        # t: translational vector
        #
        # -- out --
        #
        # v: operated u
        #
        v = [0, 0, 0]
        tt = self.tsp2frac(t)  # convert fractional notation
        for i in range(3):
            for j in range(3):
                v[i] += rot[i][j] * u[j]
            v[i] += tt[i]
        return self.reduce_coordinate(v)
