#!/usr/bin/env python3
#
# modified ParseConfig
#   written by Hiroki Funashima in Kobe 4 June 2019
#
from Metis.Base.TspaceToolbox import TspaceToolbox
import os


class ParseConfig(TspaceToolbox):
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
