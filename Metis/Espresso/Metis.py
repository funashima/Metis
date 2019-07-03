#!/usr/bin/env python3
from datetime import datetime
from Metis.Espresso.GenerateEspressoIn import GenerateEspressoIn
from Metis.SpaceGroup.GetWyckoffList import GetWyckoffList
import os
import re


class Metis(GetWyckoffList):
    def __init__(self, configfile):
        super().__init__(configfile)
        self.configfile = configfile
        self.nnode = self.configure.nnode
        self.submit_job = self.configure.submit_job
        self.main_log = 'metis.log'
        self.main()

    def generate_crystal(self, pattern):
        ispg = pattern['ispg']
        atom_info = [pattern['atom_info']]
        sub_index = pattern['sub_index']
        configfile = pattern['configfile']
        GenerateEspressoIn(configfile=configfile,
                           qe_inputfile='espresso_relax.in',
                           ispg=ispg,
                           sub_index=sub_index,
                           atom_info=atom_info,
                           submit_job=self.submit_job,
                           logfile=self.main_log)

    def main(self):
        self.write_header()
        min_spg_index = self.configure.min_spg_index
        max_spg_index = self.configure.max_spg_index
        search_pattern_lists =\
            self.setup_wyckoff_list(min_spg_index=min_spg_index,
                                    max_spg_index=max_spg_index,)

        with open(self.main_log, 'a') as fout:
            fout.write('Initially, metis found {} configurations'.
                       format(len(search_pattern_lists)))
            fout.write(' crystal structure.\n')
            fout.write('\n')
            fout.write(' metis generates crystal structure and ')
            fout.write(' check the redundancy.')
            fout.write('\n')
            fout.write('\n')
            fout.write('\n')
        for target in search_pattern_lists:
            self.generate_crystal(target)
        self.summarize()

    def write_header(self):
        with open(self.main_log, 'w') as fout:
            fout.write('==== * ==== * ==== * ====')
            fout.write(' Welcome to Metis ')
            fout.write('==== * ==== * ==== * ====')
            fout.write('\n')
            fout.write(' crystal structure search system\n')
            fout.write(' based on knowledge about ')
            fout.write(' space group and crystal symmetry\n')
            fout.write('\n')
            fout.write('           Date:{}\n'.
                       format(datetime.now().strftime("%Y/%m/%d %H:%M:%S")))
            fout.write('       Hostname:{}\n'.
                       format(os.environ['HOSTNAME']))
            fout.write('\n')
            fout.write('minimum index of space group = {}\n'.
                       format(self.configure.min_spg_index))
            fout.write('maximum index of space group = {}\n'.
                       format(self.configure.max_spg_index))
            fout.write('\n')

    def summarize(self):
        nconfig = 0
        for line in open(self.main_log, 'r'):
            linebuf = line.strip()
            if re.search('^accept', linebuf.lower()):
                nconfig += 1
        with open(self.main_log, 'a') as fout:
            print()
            print('Finally, metis found {} configurations'.
                  format(nconfig), end='')
            print(' to search crystal structure.')
            fout.write('\n')
            fout.write('Finally, metis found {} configurations'.
                       format(nconfig))
            fout.write(' to search crystal structure.')
