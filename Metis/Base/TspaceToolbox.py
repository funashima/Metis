#!/usr/bin/env python3
#
# utility and mathematical toolbox
#  written by Hiroki Funashima, 21 June 2019 in Kobe
#
# basic length unit is angstrom
#
#
from fractions import Fraction
import math


class TspaceToolbox(object):
    def get_sub_lattices(self):
        half = Fraction('1/2')
        tr1 = Fraction('1/3')
        tr2 = Fraction('1/3')
        if self.il in [0, 1]:
            pattern = []
        elif self.il == 2:
            pattern = [[0, half, half],
                       [half, 0, half],
                       [half, half, 0]]
        elif self.il == 3:
            pattern = [[half, half, half]]
        elif self.il == 4:
            pattern = [[half, half, 0]]
        elif self.il == -1:
            pattern = [[tr1, tr2, tr2],
                       [tr2, tr1, tr1]]
        else:
            print('==== error(generate_full_atomic_positions ====')
            print('il is invalid. il = {}'.format(self.il))
            exit()
        return pattern

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
        tsp_code = code.strip().lower()
        vectors = {'x': [1, 0, 0], 
                   'y': [0, 1, 0],
                   'z': [0, 0, 1],
                   'w': [1 -1, 0],
                   '-x': [-1, 0, 0], 
                   '-y': [0, -1, 0],
                   '-z': [0, 0, -1],
                   '-w': [-1, 1, 0]}
        for char, vec in vectors.items():
            if tsp_code == char:
                return vec
        self.unknown_operator

    def to_charcode(self, row):
        vectors = {'x': [1, 0, 0], 
                   'y': [0, 1, 0],
                   'z': [0, 0, 1],
                   'w': [1 -1, 0],
                   '-x': [-1, 0, 0], 
                   '-y': [0, -1, 0],
                   '-z': [0, 0, -1],
                   '-w': [-1, 1, 0]}
        for char, vec in vectors.items():
            if row == vec:
                return char
        self.unknown_operator

    def unknown_operator(self):
        print("===== error =====")
        print("unknown operator")
        exit()

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

    def ang2bohr(self, r):
        return r / 0.529177249

    def bohr2ang(self, r):
        return r * 0.529177249

    def ang2aunit(self, r):
        a = self.lattice_length[0]
        return r / a

    def bohr2aunit(self, r):
        a = self.lattice_length[0]
        return self.bohr2ang(r) / a

    def deg2rad(self, angle):
        return float(angle) * math.pi / 180.0

    def deg2cosine(self, angle):
        rad_angle = self.deg2rad(angle)
        return math.cos(rad_angle)

    def deg2sine(self, angle):
        rad_angle = self.deg2rad(angle)
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
        cn = [self.deg2cosine(x) for x in self.lattice_angle]
        sn = [self.deg2sine(x) for x in self.lattice_angle]
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

    def delta_vecs(self, u, v):
        return [v[i] - u[i] for i in range(3)]

    def init_vec(self, ndim=3):
        return [0] * ndim

    def get_nearest_dd(self, dd):
        dv = self.init_vec(ndim=len(dd))
        #
        # nearest delta_vector
        #   under considering periodic condition
        #
        for i in range(len(dd)):
            if abs(-1+dd[i]) < abs(dd[i]):
                dv[i] = -1 + dd[i]
            else:
                dv[i] = dd[i]
            if abs(1+dd[i]) < dv[i]:
                dv[i] = 1 + dd[i]
        return dv

    def _measure_length(self, u, v):
        dv = self.get_nearest_dd(self.delta_vecs(u,v))
        lvec = self.init_vec()
        for i in range(3):  # a, b, c
            for j in range(3):  # x, y, z
                lvec[j] += dv[i] * self.basis[i][j]
        return self.get_norm(lvec)

    def calc_unitcell_volume(self):
        a, b, c = self.basis
        return abs(self.inner_product(a, self.outer_product(b, c)))

    def multiply_matrices(self, r1, r2):
        #
        # calculate r1 * r2
        #
        matrix = [[0 for x in range(3)] for y in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    matrix[i][j] += r1[i][k] * r2[k][j]
        return matrix
