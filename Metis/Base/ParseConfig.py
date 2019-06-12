#!/usr/bin/env python3
#
# modified ParseConfig
#   written by Hiroki Funashima in Kobe 4 June 2019
#

import math
import os


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('file:{} is not found.'.format(configfile))
            exit()
        self.main()
        self.test_basis()

    def _get_linebuf(self, line):
        linebuf = line.strip().replace(',', ' ')
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def _get_key_and_value(self, linedata):
        key, value = [x.strip() for x in linedata.split('=')[:2]]
        key = key.lower()
        return [key, value]

    def _get_position_array(self, value):
        return [float(x) for x in self._get_linebuf(value).split()]

    def _get_value_array(self, value):
        data = [float(x) for x in value.split()]
        ndata = len(data)
        if ndata == 1:
            return [data[0]]*3
        elif ndata == 2:
            return [data[0], data[0], data[1]]
        else:
            return data

    def main(self):
        self.parse()
        self.generate_basis_vectors()

    def parse(self):
        atomic_position = False
        self.atom_list = []
        self.atomic_position = []
        self.natom = 0
        for line in open(self.configfile, 'r'):
            linebuf = self._get_linebuf(line)
            if linebuf == '':
                continue
            if atomic_position:
                if 'end_atomlist' in linebuf.lower():
                    atomic_position = False
                    continue
                data = linebuf.split(';')
                atom_name = None
                self.natom += 1
                for linedata in data:
                    if '=' in linedata:
                        key, value = self._get_key_and_value(linedata)
                    if key == 'atom_type':
                        atom_name = value
                        self.atom_list.append(atom_name)
                    elif key == 'position':
                        pos = self._get_position_array(value)
                        self.atomic_position.append(pos)
            else:
                if 'begin_atomlist' in linebuf.lower():
                    atomic_position = True
                    continue
                if '=' in linebuf:
                    key, value = self._get_key_and_value(linebuf)
                    if value == '':
                        continue
                    if key == 'lattice_length':
                        self.lattice_length = self._get_value_array(value)
                    elif key == 'lattice_angle':
                        self.lattice_angle = self._get_value_array(value)

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
        norm = 0
        for i in range(3):
            norm += v[i]*v[i]
        return math.sqrt(norm)

    def get_angle(self, u, v):
        in_prd = 0
        for i in range(3):
            in_prd += u[i] * v[i]
        # angle for cosine
        cos_ang = in_prd / (self.get_norm(u) * self.get_norm(v))
        return math.acos(cos_ang) * 180.0 / math.pi

    def _measure_length(self, u, v):
        dd = [v[i] - u[i] for i in range(3)]
        #
        # considering the periodic condition...
        #
        dv = [None, None, None]
        for i in range(3):
            dv[i] = min([abs(-1+dd[i]), abs(dd[i]), abs(1+dd[i])])
            if dd[i] < 0.0:
                dv[i] = -1.0 * dv[i]
        lvec = [0.0, 0.0, 0.0]
        for i in range(3):  # a, b, c
            for j in range(3):  # x, y, z
                lvec[j] += dv[i] * self.basis[i][j]
        return self.get_norm(lvec)

    #
    #  methods for test
    #
    def test_basis(self):
        lat = ['a', 'b', 'c']
        ang = ['alpha', 'beta', 'gamma']
        for i in range(3):
            length = self.get_norm(self.basis[i])
            err = abs(length - self.lattice_length[i])
            print('{0} = {1:8.6f} (err={2:5.3e})'.format(lat[i], length, err))

        for i in range(3):
            j1 = (i+1) % 3
            j2 = (i+2) % 3
            # u = self.basis[j1]
            # v = self.basis[j2]
            angle = self.get_angle(self.basis[j1], self.basis[j2])
            err = abs(angle - self.lattice_angle[i])
            print('{0:5s} = {1:8.5f} (err={2:5.3e})'.
                  format(ang[i], angle, err))

    def calc_length(self, i, j, aunit=False):
        #
        # self.atomic_position[k] : wyckoff position
        # self.atom_list[k] : atom_type
        #
        u = self.atomic_position[i]
        v = self.atomic_position[j]
        if aunit:
            return self._measure_length(u, v) / self.lattice_length[0]
        else:
            return self._measure_length(u, v)

