#!/usr/bin/env python3
from Metis.Espresso.QE2Spg import QE2Spg
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
import os
import re


class ParseConfig(object):
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename
        else:
            print('file:{} is not found.'.format(filename))
            exit()
        self.main()
        if os.path.isfile('espresso_relax.crystal_structure.in'):
            os.remove('espresso_relax.crystal_structure.in')

    def get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def main(self):
        self.min_spg_index = None
        self.max_spg_index = None
        self.element = None
        self.natoms = None
        for line in open(self.filename, 'r'):
            linebuf = self.get_linebuf(line)
            if linebuf == '':
                continue
            if '=' in linebuf:
                if ';' in linebuf:
                    for data in linebuf.split(';'):
                        key, value = [x.strip() for x in data.split('=')]
                        if key == 'natoms':
                            self.natoms = int(value)
                        elif key == 'element':
                            self.element = value
                else:
                    key, value = [x.strip() for x in linebuf.split('=')]
                    if value == '':
                        continue
                    if key == 'min_spg_index':
                        self.min_spg_index = int(value)
                    elif key == 'max_spg_index':
                        self.max_spg_index = int(value)


class ParseQEOutfile(object):
    def __init__(self, filename, natoms):
        if os.path.isfile(filename):
            self.filename = filename
            self.natoms = natoms
        else:
            print('file:{} is not found.'.format(filename))
            exit()
        self.main()

    def main(self):
        self.enthalpy = None
        for line in open(self.filename, 'r'):
            linebuf = line.strip()
            if 'Final' in linebuf and 'enthalpy' in linebuf:
                if re.search('^!', linebuf):
                    continue
                if '=' in linebuf:
                    self.enthalpy = float(linebuf.split('=')[1].split()[0])
                    self.enthalpy /= self.natoms
            if 'convergence' in linebuf and 'stoping' in linebuf:
                self.enthalpy = None
                return


class GetDirList(object):
    def __init__(self, min_spg_index, max_spg_index, element, natoms):
        self.min_spg_index = min_spg_index
        self.max_spg_index = max_spg_index
        self.element = element
        self.natoms = natoms
        self.main()

    def main(self):
        max_sub_index = 1000
        self.directory_list = []
        natoms_list = [x for x in range(1, self.natoms+1)
                       if self.natoms % x == 0]
        for matoms in natoms_list:
            for ispg in range(self.min_spg_index, self.max_spg_index+1):
                for sub_index in range(1, max_sub_index):
                    dirname = '{0}{1}_{2}_{3}'.format(self.element,
                                                      matoms,
                                                      ispg, sub_index)
                    if os.path.isdir(dirname):
                        self.directory_list.append([matoms, dirname])


class CheckEnthalpy(object):
    def __init__(self, configfile):
        self.configure = ParseConfig(configfile)
        self.generator = ParseGenerator()
        self.spg_dict = {}
        self.main()

    def main(self):
        self.setup()
        self.show_header()
        self.show_enthalpy()
        self.show_summary()

    def setup(self):
        self.min_spg_index = self.configure.min_spg_index
        self.max_spg_index = self.configure.max_spg_index
        self.element = self.configure.element
        self.natoms = self.configure.natoms
        self.directory_list = GetDirList(self.min_spg_index,
                                         self.max_spg_index,
                                         self.element,
                                         self.natoms).directory_list

    def show_header(self):
        print('element = {0} ; natoms = {1}'.format(self.element, self.natoms))
        print()
        print('total {} configurations'.format(len(self.directory_list)))
        print()

    def show_enthalpy(self, eps=1.0e-8):
        self.emin = None
        for (natoms, dirname) in self.directory_list:
            filename = dirname + '/espresso_relax.in'
            if not os.path.isfile(filename):
                continue
            in_spg = QE2Spg(filename).space_group.hmname
            filename = dirname + '/espresso_relax.out'
            if not os.path.isfile(filename):
                continue
            enthalpy = 0
            for line in open(filename, 'r'):
                if 'stopping' in line.strip():
                    enthalpy = None
            if enthalpy is not None:
                space_group = QE2Spg(filename).space_group.hmname
                enthalpy = ParseQEOutfile(filename, natoms).enthalpy
                #
                # CASE: enthalpy_new < enthalpy_old
                #
                if enthalpy is None:
                    for line in open(filename, 'r'):
                        if 'enthalpy old' in line.strip():
                            print(line.strip())
                            enthalpy = float(line.strip().split('=')[1].split()[0])
            if enthalpy is not None:
                if self.emin is None or self.emin > enthalpy:
                    self.emin = enthalpy
                    self.emin_spg = space_group
                    self.emin_spg_type = dirname
                if space_group in self.spg_dict:
                    check = []
                    for j in range(len(self.spg_dict[space_group])):
                        if abs(self.spg_dict[space_group][j] - enthalpy) < eps:
                            check.append(True)

                    if not all(check):
                        self.spg_dict[space_group].append(enthalpy)
                else:
                    self.spg_dict[space_group] = [enthalpy]

            print(' DIR:{0:10s} '.format(dirname), end='')
            print('spg(in):{:>6s} -> '.format(in_spg), end='')
            print('spg(out):{:>6s} ;  '.format(space_group), end='')
            if enthalpy is not None:
                print(' enthalpy H = {:12.7f} Ry per atom'.
                      format(enthalpy))
            else:
                print(' not converged case.')
                continue

    def show_summary(self):
        if self.emin is not None:
            print()
            print('---- summary ---')
            print(' we estimate most stable crystal structure,')
            print('   space_group:{0} ({1})'.
                  format(self.emin_spg, self.emin_spg_type))
            print('   min of enthalpy  = {} Ry per atom'.
                  format(self.emin))

            print()
            print('classify to each space group: delta H := H - min(H)')
            configs = 0
            for space_group, earray in self.spg_dict.items():
                print(' space_group:{:8s}  '.format(space_group), end='')
                print('delta H = ', end='')
                for enthalpy in earray:
                    configs += 1
                    print(' {:7.4e}'.format(enthalpy - self.emin), end='')
                print(' Ry per atom')
            print()
            print('actually {} configurations'.format(configs))
