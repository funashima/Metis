#!/usr/bin/env python3
#
# convert pwscf inputfile to space group
#   written by Hiroki Funashima in Kobe, 30 May 2019
#   modified by Hiroki Funashima in Kobe, 3 June 2019
#
import math


class QEFileParseBase(object):
    def __init__(self):
        pass

    def _get_linebuf(self, line):
        linebuf = line.strip()
        comment_chars = ['#', '!']
        for char in comment_chars:
            if char in linebuf:
                linebuf = linebuf.split(char)[0].strip()
        return linebuf

    def _get_key_and_value(self, keyword):
        for line in open(self.inputfile, 'r'):
            linebuf = self._get_linebuf(line)
            if linebuf == '':
                continue
            if '=' in linebuf:
                if keyword in linebuf:
                    for data in [x.strip() for x in linebuf.split(',')]:
                        if keyword in data:
                            return data.split('=')[1].strip()
        return None

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

    def _convert_cparams2lparams(self, cell_parameters):
        #
        #  convert cell_parameters -> a, b, c, alpha, beta and gamma
        #
        self.cell_length_a = self.get_norm(cell_parameters[0])
        self.cell_length_b = self.get_norm(cell_parameters[1])
        self.cell_length_c = self.get_norm(cell_parameters[2])
        self.cell_angle_alpha = self.get_angle(cell_parameters[1],
                                               cell_parameters[2])
        self.cell_angle_beta = self.get_angle(cell_parameters[2],
                                              cell_parameters[0])
        self.cell_angle_gamma = self.get_angle(cell_parameters[0],
                                               cell_parameters[1])

    def _parse_apos(self, linebuf):
        data = linebuf.split()
        atom_site_symbol = data[0]
        position = []
        for x in data[1:4]:
            position.append(eval(x))
        return [atom_site_symbol, position]

    def write_crystal_in(self):
        with open(self.outfile, 'w') as fout:
            fout.write('#\n')
            fout.write('# lattice constant\n')
            fout.write('#\n')
            fout.write('\n')
            fout.write('lattice_length = ')
            fout.write('{0:>10s}, '.format(str(self.cell_length_a)))
            fout.write('{0:>10s}, '.format(str(self.cell_length_b)))
            fout.write('{0:>10s}\n'.format(str(self.cell_length_c)))
            fout.write('lattice_angle  = ')
            fout.write('{0:>10s}, '.format(str(self.cell_angle_alpha)))
            fout.write('{0:>10s}, '.format(str(self.cell_angle_beta)))
            fout.write('{0:>10s}\n'.format(str(self.cell_angle_gamma)))
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
        return self
