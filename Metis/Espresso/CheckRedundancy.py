#!/usr/bin/env python3
import os
from Metis.SpaceGroup.GenerateSmallerUnitCell import GenerateSmallerUnitCell
import subprocess


class CheckRedundancy(object):
    #
    # return value: [consistent_spg, space_group_name]
    #

    def check(self, compound_name, ispg, sub_index):
        self.prefix = '{0}_{1}_{2}'.format(compound_name,
                                           ispg,
                                           sub_index)
        if not os.path.isdir(self.prefix):
            print('directry:{} is not found.'.format(self.prefix))
            exit()
        self.analyze_true_symmetry()
        if self.prefix == self.prim_cell.dirname:
            return [True, self.prefix]
        return [False, self.prim_cell]

    def analyze_true_symmetry(self):
        inputfile = self.prefix + '/espresso_relax.in'
        if not os.path.isfile(inputfile):
            print('===== Error(CheckRedundancy) =====')
            print('file:{} is not found.'.format(inputfile))
            subprocess.call('pwd', shell=True)
            exit()
        self.prim_cell = GenerateSmallerUnitCell(inputfile)
        if not self.prim_cell.identified_spg:
            #
            # in this case, spglib cannot identify space group for smaller prim cell
            #
            self.prim_cell.dirname = self.prefix
