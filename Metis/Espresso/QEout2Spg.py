#!/usr/bin/env python3
#
# convert pwscf inputfile to space group
#   written by Hiroki Funashima in Kobe, 30 May 2019
#   modified by Hiroki Funashima in Kobe, 3 June 2019
#
from Metis.Espresso.QEFileParseBase import QEFileParseBase
import os
import re


class QEout2Spg(QEFileParseBase):
    def __init__(self, datafile, configfile):
        if os.path.isfile(datafile):
            self.datafile = datafile
        else:
            print('file:{} is not found.'.format(datafile))
            exit()
        self.outfile = configfile
        self.main()

    def main(self):
        self.parse()
        self._convert_cparams2lparams(self.cell_parameter)
        self.length_unit_normalize()
        self.write_crystal_in()

    def length_unit_normalize(self):
        a_normalize = self.get_norm(self.cell_parameter[0])
        self.cell_length_a /= a_normalize
        self.cell_length_b /= a_normalize
        self.cell_length_c /= a_normalize

    def parse(self):
        data_region = False
        for line in open(self.datafile, 'r'):
            linebuf = self._get_linebuf(line)
            if linebuf == '':
                continue
            if data_region:
                if re.search('^density', linebuf):
                    continue
                if re.search('^new unit-cell', linebuf):
                    continue
                if re.search('end final coordinates', linebuf.lower()):
                    data_region = False
                    continue
                if re.search('CELL_PARAMETER', linebuf):
                    data_type = 'cell_parameter'
                    self.alat = float(linebuf.replace(')', '').split('=')[-1])
                    continue
                elif re.search('ATOMIC_POSITIONS', linebuf):
                    data_type = 'atomic_positions'
                    continue
                if data_type == 'cell_parameter':
                    cell_params = [float(x) for x in linebuf.split()]
                    self.cell_parameter.append(cell_params)
                elif data_type == 'atomic_positions':
                    atom_site_symbol, position = self._parse_apos(linebuf)
                    info = {'site_type': atom_site_symbol,
                            'coordinate': position}
                    self.atom_site.append(info)
            else:
                if re.search('begin final coordinates', linebuf.lower()):
                    data_region = True
                    self.cell_parameter = []
                    self.atom_site = []
                    continue
