#!/usr/bin/env python3
import os
import re


class ParseTestInput(object):
    def __init__(self, inputfile):
        if os.path.isfile(inputfile):
            self.inputfile = inputfile
        else:
            print('file:{} is not found.'.format(inputfile))
            exit()
        self.main()

    def main(self):
        atom_info_region = False
        self.ichoice = 1
        self.max_coa_ratio = 2.0
        self.atom_info = []
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
                wypos = None
                for data in linebuf.split(';')[:2]:
                    key, value = [x.strip() for x in data.split('=')[:2]]
                    key = key.lower()
                    if key == 'element':
                        element = value
                    if key == 'wyckoff_position':
                        wypos = value.replace(',', ' ').split()
                if (element is not None) and (wypos is not None):
                    info = {'element': element, 'wyckoff_position': wypos}
                    self.atom_info.append(info)
            else:
                if re.search('^begin_atom_info:', linebuf):
                    if re.search('^begin_atom_info:', linebuf):
                        atom_info_region = True
                        continue
                if '=' in linebuf:
                    key, value = [x.strip() for x in linebuf.split('=')[:2]]
                    key = key.lower()
                    if key == 'space_group':
                        self.space_group = value
                    if key == 'ichoice':
                        self.ichoice = int(value)
                    if key == 'max_coa_ratio':
                        self.max_coa_ratio = float(value)
                    if key == 'apf':
                        self.apf = float(value)
                    if key == 'max_try':
                        self.max_try = int(value)
