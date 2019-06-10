#!/usr/bin/env python3
#
# cif2space group
#   written by Hiroki Funashima in Kobe, 29 May 2019
#

import os
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import re


class Cif2Spg(object):
    def __init__(self, ciffile='crystal.cif', outfile='crystal.in'):
        if os.path.isfile(ciffile):
            self.ciffile = ciffile
        else:
            print('===== Error(cif2spg) ====')
            print('file:{} is not found'.format(ciffile))
            exit()
        self.outfile = outfile
        self.parse()
        self.write_crystal_in()

    def _get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def _get_key_and_value(self, line):
        linebuf = self._get_linebuf(line)
        if not re.search('^_', linebuf):
            return [None, None]
        data = linebuf.split()
        if len(data) < 2:
            key = data[0]
            value = None
        else:
            key, value = linebuf.split()[:2]
            if value == '':
                value = None
            if '(' in value:
                value = value.split('(')[0].strip()
            if re.search('.$', value):
                value += '0'
        return [key, value]

    def my_split(self, s, sep):
        kakkoDic = {'"': '"',
                    '\'': '\'',
                    '[': ']'}
        kakkoKeys = kakkoDic.keys()
        reslist = []
        ary = s.split(sep)
        targetstr = ary[0]
        if len(targetstr) >= 1 and targetstr[0] in kakkoKeys:
            p = re.compile(kakkoDic[targetstr[0]]+'$')
            joinlist = []
            num = 0
            for s in ary:
                num = num + 1
                joinlist.append(s)
                if len(ary) <= 1 or p.search(s) is not None:
                    s2 = joinlist[len(joinlist) - 1]
                    joinlist[len(joinlist)-1] = s2[0:len(s2)]
                    break
            c = sep.join(joinlist)
            c = c[1:len(c) - 1]
            ary = ary[num:len(ary)]
            reslist.append(c)
        else:
            reslist.append(ary[0])
            ary.pop(0)
        if ary:
            reslist = reslist + self.my_split(sep.join(ary), sep)
        return [x for x in reslist if x.strip() != '']

    def parse(self):
        self.get_lattice_param()
        self.get_symmetry_operation()
        self.get_atomic_position()
        self.apply_symmetry_operation()

    def initialize_lattice_params(self):
        self.lattice_params = ['cell_length_a',
                               'cell_length_b',
                               'cell_length_c',
                               'cell_angle_alpha',
                               'cell_angle_beta',
                               'cell_angle_gamma']
        for variable in self.lattice_params:
            exec('self.{} = None'.format(variable))

    def get_lattice_param(self):
        self.initialize_lattice_params()
        for line in open(self.ciffile, 'r'):
            linebuf = self._get_linebuf(line)
            if linebuf == '':
                continue
            key, value = self._get_key_and_value(line)
            if value is not None:
                for variable in self.lattice_params:
                    keyword = '_{}'.format(variable)
                    if key == keyword:
                        exec('self.{0} = "{1}"'.format(variable, value))

    def get_loop_index(self, keyword):
        loop_db = []
        loop_check = False
        for line in open(self.ciffile, 'r'):
            linebuf = self._get_linebuf(line)
            if 'loop_' in linebuf:
                loop_check = True
                loop_db.append([])
            else:
                if loop_check:
                    if re.search('^_', linebuf):
                        title = linebuf.split()[0]
                        loop_db[-1].append(title)
        for list_name in loop_db:
            if keyword in list_name:
                return list_name.index(keyword)

    def get_symmetry_operation(self):
        data_region = False
        self.symmetry_operation = []
        keyword = '_symmetry_equiv_pos_as_xyz'
        frac_id = self.get_loop_index(keyword)

        for line in open(self.ciffile):
            linebuf = self._get_linebuf(line)
            if data_region:
                if linebuf == '' or linebuf == 'loop_':
                    data_region = False
                    continue
                if re.search('^_', linebuf):
                    continue
                word = self.my_split(linebuf, " ")[frac_id]
                data = [x.strip() for x in word.replace("'", '').split(",")]
                self.symmetry_operation.append(data)
            else:
                if keyword in linebuf:
                    data_region = True

    def get_atomic_position(self):
        self.atomic_position = []
        data_region = False
        x_id = self.get_loop_index('_atom_site_fract_x')
        y_id = self.get_loop_index('_atom_site_fract_y')
        z_id = self.get_loop_index('_atom_site_fract_z')
        label_id = self.get_loop_index('_atom_site_type_symbol')
        for line in open(self.ciffile, 'r'):
            linebuf = self._get_linebuf(line)
            if data_region:
                if 'loop_' in linebuf or linebuf == '':
                    data_region = False
                    continue
                if not re.search('^_', linebuf):
                    data = linebuf.split()
                    atom_site_label = data[label_id]
                    atom_site_fract_x = data[x_id]
                    if '(' in atom_site_fract_x:
                        atom_site_fract_x = atom_site_fract_x.split('(')[0]
                    atom_site_fract_y = data[y_id]
                    if '(' in atom_site_fract_y:
                        atom_site_fract_y = atom_site_fract_y.split('(')[0]
                    atom_site_fract_z = data[z_id]
                    if '(' in atom_site_fract_z:
                        atom_site_fract_z = atom_site_fract_z.split('(')[0]
                    coordinate = [atom_site_fract_x,
                                  atom_site_fract_y,
                                  atom_site_fract_z]
                    pos_data = {'atom_site_type_symbol': atom_site_label,
                                'coordinate': coordinate}
                    self.atomic_position.append(pos_data)
            if '_atom_site_fract_x' in linebuf:
                data_region = True
                continue

    def _renew_coordinate(self, position):
        pos = [None] * 3
        for i in range(3):
            pos[i] = self._reduce_bz(position[i])
        return pos

    def _reduce_bz(self, x):
        while abs(x) > 0.5:
            if x > 0.5:
                x -= 1.0
            if x < -0.5:
                x += 1.0
        return x

    def equivalent_r(self, label0, r0,  label1, r1, eps=1.0e-6):
        if label0 != label1:
            return False
        for i in range(3):
            if abs(r0[i] - r1[i]) > eps:
                return False
        return True

    def apply_symmetry_operation(self):
        self.atom_site = []
        for atom_data in self.atomic_position:
            position = atom_data['coordinate']
            xx = [None] * 3
            for j in range(3):
                xx[j] = '({})'.format(position[j])
            label = atom_data['atom_site_type_symbol']
            for sym_op in self.symmetry_operation:
                atomic_position = [None, None, None]
                for i in range(3):
                    position = sym_op[i]
                    x = ['x', 'y', 'z']
                    for j in range(3):
                        position = position.replace(x[j], xx[j])
                    try:
                        atomic_position[i] = eval(position)
                    except SyntaxError:
                        print('eval error in apply_symmetry_operation')
                        print(position)
                        exit()
                atomic_position = self._renew_coordinate(atomic_position)

                label1 = label
                r1 = atomic_position
                info = {'site_type': label1, 'coordinate': r1}
                if len(self.atom_site) < 1:
                    self.atom_site.append(info)
                else:
                    reg_check = True
                    for registered in self.atom_site:
                        label0 = registered['site_type']
                        r0 = registered['coordinate']
                        if self.equivalent_r(label0, r0,  label1, r1):
                            reg_check = False
                    if reg_check:
                        self.atom_site.append(info)

    def show_info(self):
        print('-- lattice params --')
        print(' a = {}'.format(self.cell_length_a))
        print(' b = {}'.format(self.cell_length_b))
        print(' c = {}'.format(self.cell_length_c))
        print(' alpha = {}'.format(self.cell_angle_alpha))
        print(' beta  = {}'.format(self.cell_angle_beta))
        print(' gamma = {}'.format(self.cell_angle_gamma))
        print()
        print(self.symmetry_operation)
        print(self.atomic_position)

    def write_crystal_in(self):
        with open(self.outfile, 'w') as fout:
            fout.write('#\n')
            fout.write('# lattice constant\n')
            fout.write('#\n')
            fout.write('\n')
            fout.write('lattice_length = ')
            fout.write('{0:>10s}, '.format(self.cell_length_a))
            fout.write('{0:>10s}, '.format(self.cell_length_b))
            fout.write('{0:>10s}\n'.format(self.cell_length_c))
            fout.write('lattice_angle  = ')
            fout.write('{0:>10s}, '.format(self.cell_angle_alpha))
            fout.write('{0:>10s}, '.format(self.cell_angle_beta))
            fout.write('{0:>10s}\n'.format(self.cell_angle_gamma))
            fout.write('\n')
            fout.write('#\n')
            fout.write('# atomic position in conventional unit cell\n')
            fout.write('#\n')
            fout.write('begin_atomlist:\n')
            for info in self.atom_site:
                label = info['site_type']
                atomic_position = info['coordinate']
                fout.write('    ')
                fout.write('atom_type = {} ; position = '.format(label))
                for i in range(3):
                    fout.write('{0:10.7f}'.format(atomic_position[i]))
                    if i < 2:
                        fout.write(', ')
                fout.write('\n')
            fout.write('end_atomlist:\n')


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('file:{} is not found.'.format(configfile))
            exit()
        self.main()

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
        atomic_position = False
        self.atom_list = []
        self.atomic_position = []
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


class SpaceGroup(object):
    def __init__(self, configfile, symprec=0.01,  angle_trelance=6):
        self.config_obj = ParseConfig(configfile)
        self.symprec = symprec
        self.angle_trelance = angle_trelance
        self.main()

    def set_structure(self):
        #
        # ======= ATOMIC POSITIONS =======
        #
        a, b, c = self.config_obj.lattice_length[:3]
        alpha, beta, gamma = self.config_obj.lattice_angle[:3]
        lattice = mg.Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        self.atoms = self.config_obj.atom_list
        atomic_positions = self.config_obj.atomic_position
        self.structure = mg.Structure(lattice, self.atoms, atomic_positions)

    def main(self):
        self.set_structure()
        #
        # generate instance for space group
        #
        self.spg = SpacegroupAnalyzer(self.structure,
                                      symprec=self.symprec,
                                      angle_tolerance=self.angle_trelance)

    def show_info(self):
        print('----- symmetrized structure(conventional unit cell) -----')
        print(self.spg.get_symmetrized_structure())
        print()
        print('----- primitive structure -----')
        print(self.spg.find_primitive())
        print()
        #
        # name of space group
        #
        print('--- infomation of space group ---')
        print('  HM symbol: {} '.format(self.spg.get_space_group_symbol()))
        print('  Space Group number = #{}'.
              format(self.spg.get_space_group_number()))
        print('  point group : {}'.format(self.spg.get_point_group_symbol()))
        print()
        print()
        print('--- crystal system ---')
        print('  crystal system:{}'.format(self.spg.get_crystal_system()))
        print('  lattice type:{}'.format(self.spg.get_lattice_type()))

        #
        # total data set
        #
        dataset = self.spg.get_symmetry_dataset()
        wyckoff_position = dataset['wyckoffs']
        atom_pos = []
        atom_name = []
        wyckoff_data = []
        for (i, wyckoff) in enumerate(wyckoff_position):
            aname = self.atoms[i]
            if aname not in atom_name or wyckoff not in atom_pos:
                info = {'site_type_symbol': aname,
                        'wyckoff_letter': wyckoff,
                        'multiply': 0}
                wyckoff_data.append(info)
                atom_pos.append(wyckoff)
                atom_name.append(aname)
        for (i, wyckoff) in enumerate(wyckoff_position):
            aname = self.atoms[i]
            for info in wyckoff_data:
                if info['site_type_symbol'] == aname:
                    if info['wyckoff_letter'] == wyckoff:
                        info['multiply'] += 1
        print('  # of kinds of atoms = {}'.format(len(wyckoff_data)))
        print()
        for info in wyckoff_data:
            site_name = '{0}{1}-site'.format(info['multiply'],
                                             info['wyckoff_letter'])
            atom_name = info['site_type_symbol']
            print('  atom:{0}  wyckoff: {1}'.format(atom_name, site_name))
        return

    def symmetrized(self):
        self.get_symmetrized_structure()
