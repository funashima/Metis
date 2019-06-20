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
        self.generate_full_atomic_positions()
        self.set_lattice()

    def calc_sphere_volume(self, r):
        return (4.0 * math.pi / 3.0) * (r**3)

    def get_radius(self, atomic_symbol):
        return float(mg.Element(atomic_symbol).atomic_radius)

    def set_const_volume(self, afp=0.52):
        #
        # -Atomic Packing Factor(APF)-
        #
        # APF: bcc = pi*sqrt(3.0) / 8 = 0.68
        #      fcc/hcp = pi/sqrt(18) = 0.74
        #      diamond = 0.34
        #      simple cubic lattice = 0.52
        #
        volume = 0.0
        for atom in self.atom_data:
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
        volume /= afp
        return volume

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
        self.atom_data = []
        for atom in self.atom_info:
            atom_name = atom['element']
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
            wyckoff_letter = atom_site.wyckoff_params['letter']
            self.atom_data.append({'element': atom_name_list[i],
                                   'wyckoff_position': wyckoff_letter,
                                   'positions': atomic_pos})

    def ang2bohr(self, r):
        return r / 0.529177249

    def ang2aunit(self, r):
        a = self.config_obj.lattice_length[0]
        return self.ang2bohr(r) / a

    def _deg2rad(self, angle):
        return float(angle) * math.pi / 180.0

    def _deg2cosine(self, angle):
        rad_angle = self._deg2rad(angle)
        return math.cos(rad_angle)

    def _deg2sine(self, angle):
        rad_angle = self._deg2rad(angle)
        return math.sin(rad_angle)

    def _get_trigonometric_funcs(self):
        #
        # cn[0] = cos(alpha)
        # cn[1] = cos(beta)
        # cn[2] = cos(gamma)
        #
        # sn[0] = sin(alpha)
        # sn[1] = sin(beta)
        # sn[2] = sin(gamma)
        #
        cn = [self._deg2cosine(x) for x in self.lattice_angle]
        sn = [self._deg2sine(x) for x in self.lattice_angle]
        return [cn, sn]

    def set_basis(self, a, b, c):
        basis = [None] * 3
        cn, sn = self._get_trigonometric_funcs()
        basis[0] = [float(a), 0.0, 0.0]
        basis[1] = [a*cn[2], b*sn[2], 0.0]

        cx = cn[1]
        cy = (1/sn[2])*(cn[0] - cn[1]*cn[2])
        cz = math.sqrt(1.0 - (cx**2 + cy**2))
        basis[2] = [a * cx,  b * cy, c * cz]
        return basis

    def generate_basis_vectors(self):
        #
        a, b, c = [float(x) for x in self.lattice_length[0:3]]
        self.basis = self.set_basis(a, b, c)

    def get_norm(self, v):
        return math.sqrt(self.inner_product(v, v))

    def get_angle(self, u, v):
        in_prd = 0
        for i in range(3):
            in_prd += u[i] * v[i]
        # angle for cosine
        cos_ang = in_prd / (self.get_norm(u) * self.get_norm(v))
        return math.acos(cos_ang) * 180.0 / math.pi

    def inner_product(self, u, v):
        iprd = 0
        for i in range(3):
            iprd += u[i] * v[i]
        return iprd

    def _outer_prd_index(self, i):
        return[(i + 1) % 3, (i + 2) % 3]

    def outer_product(self, u, v):
        w = [0, 0, 0]
        for i in range(3):
            j, k = self._outer_prd_index(i)
            w[i] = u[j] * v[k] - u[k] * v[j]
        return w

    def _measure_length(self, u, v):
        dd = [v[i] - u[i] for i in range(3)]
        #
        # considering the periodic condition...
        #
        dv = [None, None, None]
        for i in range(3):
            if abs(-1+dd[i]) < abs(dd[i]):
                dv[i] = -1 + dd[i]
            else:
                dv[i] = dd[i]
            if abs(1+dd[i]) < dv[i]:
                dv[i] = 1 + dd[i]
        lvec = [0.0, 0.0, 0.0]
        for i in range(3):  # a, b, c
            for j in range(3):  # x, y, z
                lvec[j] += dv[i] * self.basis[i][j]
        return self.get_norm(lvec)

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
        for atom in self.atom_data:
            element = atom['element']
            wyckoff_position = atom['wyckoff_position']
            print('Atom:{}  wyckoff_position:{}'.
                  format(element, wyckoff_position))
            for atomic_position in atom['positions']:
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
        print('Schoenflies symbol    :{}'.
              format(self.space_group.schname))
        print('Hermann-Mauguin symbol:{}'.
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
                print('  to atom({0:>2d}): {1:>2s}'.format(j, atom_name1), end='')
                print(' length = {0:8.5f} ang'.format(length0), end='')
                print(' (= {0:6.5f})'.format(length1), end='')
                # sum of atomic radius
                r_atoms = self.get_radius(atom_name0) + self.get_radius(atom_name1)
                print('; sum of 2 atomic radii = {:7.5f} ang'.format(r_atoms), end='')
                print('  state ... ', end='')
                diff = length0 - r_atoms
                ratio = (length0 / r_atoms) * 100
                if r_atoms >= length0:
                    print('NG', end='')
                else:
                    print('OK', end='')
                print(' ({0:8.5f} ang, {1:6.2f}%)'.format(diff, ratio))


    def calc_unitcell_volume(self):
        a, b, c = self.basis
        return abs(self.inner_product(a, self.outer_product(b, c)))

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

    def add_sublattice(self, u, sub_lattice):
        v = [0, 0, 0]
        for i in range(3):
            v[i] = u[i] + sub_lattice[i]
        return self.reduce_coordinate(v)

    def generate_full_atomic_positions(self):
        half = Fraction('1/2')
        tr1 = Fraction('1/3')
        tr2 = Fraction('1/3')
        self.full_atomic_position = []
        if self.il in [0, 1]:
            self.sub_lattice_pattern = []
        elif self.il == 2:
            self.sub_lattice_pattern = [[0, half, half],
                                        [half, 0, half],
                                        [half, half, 0]]
        elif self.il == 3:
            self.sub_lattice_pattern = [[half, half, half]]
        elif self.il == 4:
            self.sub_lattice_pattern = [[half, half, 0]]
        elif self.il == -1:
            self.sub_lattice_pattern = [[tr1, tr2, tr2],
                                        [tr2, tr1, tr1]]
        else:
            print('==== error(generate_full_atomic_positions ====')
            print('il is invalid. il = {}'.format(self.il))
            exit()

        for atom in self.atom_data:
            element = atom['element']
            wyckoff_position = atom['wyckoff_position']
            for v in atom['positions']:
                info = {'element': element,
                        'wyckoff_position': wyckoff_position,
                        'positions':v}
                self.full_atomic_position.append(info)
                for sub_lattice in self.sub_lattice_pattern:
                    vv = self.add_sublattice(v, sub_lattice)
                    info = {'element': element,
                            'wyckoff_position': wyckoff_position,
                            'positions': vv}
                    self.full_atomic_position.append(info)

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
