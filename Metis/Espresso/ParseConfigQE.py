#!/usr/bin/env python3
#
# 
#
import os
import re


class ParseConfigQE(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('===== Error(ParseConfig) =====')
            print('file:{} is not found.'.format(configfile))
            exit()
        self.main()

    def get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def get_key_and_value(self, linebuf):
        if '=' not in linebuf:
            return [None, None]
        key, value = [x.strip() for x in linebuf.split('=')[:2]]
        key = key.lower()
        return [key, value]

    def value2bool(self, value):
        if re.search('^(t|y)', value.lower()):
            return True
        elif re.search('^(f|n)', value.lower()):
            return False
        return False

    def check_pslib(self, value):
        if os.path.isdir(value):
            self.pslib = value
        else:
            print('===== Error(ParseConfig) =====')
            print('dir:{} is not found'.format(value))
            exit()

    def set_environments(self, linebuf):
        key, value = self.get_key_and_value(linebuf)
        if key is None:
            return
        if value == '':
            return
        if key == 'pslib':
            self.check_pslib(value)
        elif key == 'dft_type':
            self.dft_type = value
        elif key == 'pp_type':
            self.pp_type = value
        elif key == 'spin_orbit':
            self.spin_orbit = self.value2bool(value)
        elif key == 'degauss':
            self.degauss = float(value)
        elif key == 'electron_maxstep':
            self.electron_maxstep = int(value)
        elif key == 'conv_thr':
            self.conv_thr = float(value)
        elif key == 'mixing_beta':
            self.mixing_beta = float(value)
        elif key == 'kpoints':
            self.kpoints = [int(x) for x in value.replace(',', ' ').split()]
            if len(self.kpoints) == 1:
                kp = self.kpoints[0]
                self.kpoints = [kp, kp, kp]
            elif len(self.kpoints) == 2:
                kp1 = self.kpoints[0]
                kp3 = self.kpoints[1]
                self.kpoints = [kp1, kp1, kp3]
        elif key == 'ecutwfc':
            self.ecutwfc = float(value)
        elif key == 'ecutrho':
            self.ecutrho = float(value)
        elif key == 'use_auto_cutoff':
            if re.search('(t|y)', value.lower()):
                self.use_auto_cutoff = True
            else:
                self.use_auto_cutoff = False


        elif key == 'press':
            self.press = float(value)
        elif key == 'cell_factor':
            self.cell_factor = float(value)

    def set_atom_info(self, linebuf):
        info = {'element': None,
                'natoms': 1,
                'semi_core': []}
        for linedata in [x.strip() for x in linebuf.split(';')]:
            key, value = self.get_key_and_value(linedata)
            if key is None:
                continue
            if key == 'element':
                info['element'] = value
            if key == 'natoms':
                info['natoms'] = int(value)
            if key == 'semi_core':
                info['semi_core'] = value.replace(',', ' ').split()
        if info['element'] is not None:
            self.atom_info.append(info)

    def set_init_value(self):
        self.dft_type = 'pbe'
        self.pp_type = 'USPP'
        self.spin_orbit = False
        self.atom_info = []
        self.degauss = 0.02
        self.electron_maxstep = 100
        self.conv_thr = 1.0e-6
        self.mixing_beta = 0.7
        self.kpoints = [6, 6, 6]
        self.press = 0.0
        self.cell_factor = 1.0
        self.ecutwfc = None
        self.ecutrho = None
        self.use_auto_cutoff = False

    def main(self):
        self.set_init_value()
        atom_info_region = False
        for line in open(self.configfile, 'r'):
            linebuf = self.get_linebuf(line)
            if linebuf == '':
                continue
            if atom_info_region:
                if re.search('^end_atom_info', linebuf.lower()):
                    atom_info_region = False
                    continue
                self.set_atom_info(linebuf)
            else:
                if re.search('^begin_atom_info', linebuf.lower()):
                    atom_info_region = True
                    continue
                if '=' in linebuf:
                    self.set_environments(linebuf)
