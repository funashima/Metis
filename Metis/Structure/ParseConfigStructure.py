#!/usr/bin/env python3
import os
import re


class ParseConfigStructure(object):
    def __init__(self, inputfile):
        if os.path.isfile(inputfile):
            self.inputfile = inputfile
        else:
            print('file:{} is not found.'.format(inputfile))
            exit()
        self.main()

    def set_default_value(self):
        self.min_spg_index = 8
        self.max_spg_index = 230
        self.ichoice = 1
        self.max_coa_ratio = 2.0
        self.atom_info = []
        self.apf = 1.0
        self.delata_apf = 0.9
        self.thr_bond_ratio = 0.75
        self.max_try = 500
        self.nnode = 1

    def main(self):
        atom_info_region = False
        self.set_default_value()
        for line in open(self.inputfile, 'r'):
            linebuf = line.strip()
            if '#' in linebuf:
                linebuf = linebuf.split('#')[0].strip()
            if linebuf == '':
                continue
            if atom_info_region:
                if re.search('^end_atom_info:', linebuf):
                    atom_info_region = False
                    continue
                element = None
                natoms = None
                semi_core = []
                for data in linebuf.split(';'):
                    if '=' in data:
                        if len(data.split('=')) > 1:
                            key, value = [x.strip() for x in data.split('=')[:2]]
                            key = key.lower()
                            if key == 'element':
                                element = value
                            if key == 'natoms':
                                natoms = int(value)
                            if key == 'semi_core':
                                semi_core = value.replace(',', ' ').split()
                if (element is not None) and (natoms is not None):
                    info = {'element': element, 'natoms': natoms, 'semi_core': semi_core}
                    self.atom_info.append(info)
                else:
                    print('WARNING: syntax is incorrect :{}'.format(linebuf))
            else:
                if re.search('^begin_atom_info:', linebuf):
                    if re.search('^begin_atom_info:', linebuf):
                        atom_info_region = True
                        continue
                if '=' in linebuf:
                    if len(linebuf.split('=')) < 2:
                        continue
                    key, value = [x.strip() for x in linebuf.split('=')[:2]]
                    key = key.lower()
                    if key == 'min_spg_index':
                        self.min_spg_index = int(value)
                    if key == 'max_spg_index':
                        self.max_spg_index = int(value)
                    elif key == 'ichoice':
                        self.ichoice = int(value)
                    elif key == 'max_coa_ratio':
                        self.max_coa_ratio = float(value)
                    elif key == 'apf':
                        self.apf = float(value)
                    elif key == 'max_try':
                        self.max_try = int(value)
                    elif key == 'delta_apf':
                        self.delta_apf = float(value)
                    elif key == 'thr_bond_ratio':
                        self.thr_bond_ratio = float(value)
                    elif key == 'nnode':
                        self.nnode = int(value)
