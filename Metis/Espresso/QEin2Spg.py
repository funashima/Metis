#!/usr/bin/env python3
#
# convert pwscf inputfile to space group
#   written by Hiroki Funashima in Kobe, 30 May 2019
#   modified by Hiroki Funashima in Kobe, 3 June 2019
#

from Metis.Espresso.QEFileParseBase import QEFileParseBase
import math
import os


class QEin2Spg(QEFileParseBase):
    def __init__(self, inputfile, outfile):
        if os.path.isfile(inputfile):
            self.inputfile = inputfile
        else:
            print('===== Error(QEin2Spg) =====')
            print('file:{} is not found.'.format(inputfile))
            exit()
        self.outfile = outfile
        self.parse()
        self.write_crystal_in()

    def parse(self):
        self.get_lattice_params()
        self.get_cell_parameters()
        self.get_atomic_positions()

    def get_nat(self):
        self.nat = self._get_key_and_value('nat')
        if self.nat is not None:
            self.nat = int(self.nat)

    def get_ibrav(self):
        self.ibrav = self._get_key_and_value('ibrav')
        if self.ibrav is not None:
            self.ibrav = int(self.ibrav)

    def get_cell_parameters(self):
        cell_parameters = []
        data_region = False
        ix = 0
        for line in open(self.inputfile, 'r'):
            linebuf = self._get_linebuf(line)
            if data_region:
                ix += 1
                if ix == 4:
                    data_region = False
                    continue
                x = [float(x) for x in linebuf.split()[:3]]
                cell_parameters.append(x)
            else:
                if 'CELL_PARAMETERS' in linebuf.upper():
                    data_region = True
                    continue
        if len(cell_parameters) == 3:
            self._convert_cparams2lparams(cell_parameters)
        return self

    def init_celldm(self):
        for i in range(1, 7):
            exec('self.celldm{} = None'.format(i))
        return self

    def set_cell_length_and_angle(self):
        self.cell_length_a = self.celldm1
        if self.ibrav == 1 or self.ibrav == 2 or self.ibrav == 3:
            self.cell_length_b = self.cell_length_a
            self.cell_length_c = self.cell_length_a
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = 90.0
            self.cell_angle_gamma = 90.0
        elif self.ibrav == 4:  # hexagonal and rhombohedral(hex)
            self.cell_length_b = self.cell_length_a
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = 90.0
            self.cell_angle_gamma = 120.0
        elif self.ibrav == 5:  # rhombohedral
            self.cell_length_b = self.cell_length_a
            self.cell_length_c = self.cell_length_a
            self.cell_angle_alpha = math.acos(self.celldm4) * 180.0 / math.pi
            self.cell_angle_beta = self.cell_angle_alpha
            self.cell_angle_gamma = self.cell_angle_alpha
        elif self.ibrav == 6 and self.ibrav == 7:  # tetragonal
            self.cell_length_b = self.cell_length_a
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = 90.0
            self.cell_angle_gamma = 90.0
        elif 8 <= self.ibrav <= 11:  # orthorhombic
            self.cell_length_b = self.cell_length_a * self.celldm2
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = 90.0
            self.cell_angle_gamma = 90.0
        elif self.ibrav == 12 and self.ibrav == 13:  # monoclinic
            self.cell_length_b = self.cell_length_a * self.celldm2
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = 90.0
            self.cell_angle_gamma = math.acos(self.celldm4) * 180.0 / math.pi
        elif self.ibrav == -12 and self.ibrav == -13:  # monoclinic
            self.cell_length_b = self.cell_length_a * self.celldm2
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = 90.0
            self.cell_angle_beta = math.acos(self.celldm5) * 180.0 / math.pi
            self.cell_angle_gamma = 90.0
        elif self.ibrav == 14:  # triclinic
            self.cell_length_b = self.cell_length_a * self.celldm2
            self.cell_length_c = self.cell_length_a * self.celldm3
            self.cell_angle_alpha = math.acos(self.celldm4) * 180.0 / math.pi
            self.cell_angle_beta = math.acos(self.celldm5) * 180.0 / math.pi
            self.cell_angle_gamma = math.acos(self.celldm6) * 180.0 / math.pi
        else:
            print('not implement type ibrav = {}'.format(self.ibrav))
            exit()

        self.set_il()
        return self

    def set_il(self):
        if self.ibrav in [0, 1, 4, 5, 6, 8, 12, -12, 14]:
            self.il = 'P'
        if self.ibrav == 2:
            self.il = 'F'
        if self.ibrav in [3, -3, 7, 11]:
            self.il = 'I'
        if self.ibrav in [9]:
            self.il = 'C'
        return self

    def get_lattice_params(self):
        self.get_ibrav()
        self.init_celldm()
        if self.ibrav == 0:
            self.il = 'P'
            self.get_cell_parameters()
        for i in range(1, 7):
            keyword = 'celldm({})'.format(i)
            value = self._get_key_and_value(keyword)
            if value is not None:
                try:
                    exec('self.celldm{0} = {1}'.format(i, value))
                except SyntaxError:
                    print('self.celldm{0} = {1}'.format(i, value))
                    exit()
        if self.ibrav != 0:
            if self.celldm1 is None:
                self.get_lattice_param_type2()
            else:
                self.set_cell_length_and_angle()

    def get_lattice_param_type2(self):
        self.cell_length_a = self._get_key_and_value('A')
        self.cell_length_b = self.cell_length_a
        self.cell_length_c = self.cell_length_a
        self.cell_angle_alpha = 90.0
        self.cell_angle_beta = 90.0
        if self.ibrav == 4:
            self.cell_angle_gamma = 120.0
        else:
            self.cell_angle_gamma = 90

        self.cell_length_b = self._get_key_and_value('B')
        self.cell_length_c = self._get_key_and_value('C')

        value = self._get_key_and_value('cosBC')
        if value is not None:
            self.cell_angle_alpha = math.acos(float(value))

        value = self._get_key_and_value('cosAC')
        if value is not None:
            self.cell_angle_beta = math.acos(float(value))

        value = self._get_key_and_value('cosAB')
        if value is not None:
            self.cell_angle_gamma = math.acos(float(value))
        self.set_il()

    def _set_atomic_position_coefficient(self, linebuf):
        xtype = None
        if '{' in linebuf:
            linebuf = linebuf.replace('{', ' ').replace('}', '')
        data = linebuf.split()
        if len(data) > 1:
            xtype = data[1].lower()
        else:
            xtype = 'alat'
        self.apos = self._xtype2coefficient(xtype)

    def _xtype2coefficient(self, xtype):
        if xtype == 'alat':
            return [1.0,
                    (self.cell_length_a / self.cell_length_b),
                    (self.cell_length_a / self.cell_length_c)]
        elif xtype == 'crystal':
            return [1.0, 1.0, 1.0]
        elif xtype == 'angstrom':
            return [1.0 / (self.cell_length_a * 0.529177249),
                    1.0 / (self.cell_length_b * 0.529177249),
                    1.0 / (self.cell_length_c * 0.529177249)]
        elif xtype == 'bohr':
            return [1.0 / (self.cell_length_a),
                    1.0 / (self.cell_length_b),
                    1.0 / (self.cell_length_c)]
        else:
            print('not implement type:{}'.format(xtype))
            exit()

    def get_atomic_positions(self):
        data_region = False
        self.get_nat()
        iatom = 0
        self.atom_site = []
        for line in open(self.inputfile, 'r'):
            linebuf = self._get_linebuf(line)
            if linebuf == '':
                continue
            if data_region:
                iatom += 1
                if iatom > self.nat:
                    data_region = False
                    continue
                atom_site_symbol, position = self._parse_apos(linebuf)
                info = {'site_type': atom_site_symbol,
                        'coordinate': position}
                self.atom_site.append(info)
                if self.il == 'P':
                    continue
                if self.il == 'F':
                    sub = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
                elif self.il == 'I':
                    sub = [[0.5, 0.5, 0.5]]
                elif self.il == 'C':
                    sub = [[0.5, 0.5, 0.0]]
                else:
                    print('il error')
                    exit()
                for lat in sub:
                    pos = [None, None, None]
                    for i in range(3):
                        pos[i] = position[i] + lat[i]
                    info = {'site_type': atom_site_symbol,
                            'coordinate': pos}
                    self.atom_site.append(info)

            else:
                if 'ATOMIC_POSITIONS' in linebuf.upper():
                    self._set_atomic_position_coefficient(linebuf)
                    data_region = True
                    continue

    def show_info(self):
        print('natom = {}'.format(self.nat))
        print('ibrav = {}'.format(self.ibrav))
        print('il = {}'.format(self.il))
        print('cell_length_a = {}'.format(self.cell_length_a))
        print('cell_length_b = {}'.format(self.cell_length_b))
        print('cell_length_c = {}'.format(self.cell_length_c))
        print('cell_angle_alpha = {}'.format(self.cell_angle_alpha))
        print('cell_angle_beta  = {}'.format(self.cell_angle_beta))
        print('cell_angle_gamma = {}'.format(self.cell_angle_gamma))
        for i in range(1, 7):
            celldm = eval('self.celldm{}'.format(i))
            print('celldm({0}) = {1}'.format(i, celldm))
        for info in self.atom_site:
            print(info)
