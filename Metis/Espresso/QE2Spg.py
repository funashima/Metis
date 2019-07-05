#!/usr/bin/env python3
#
# convert pwscf inputfile to space group
#   written by Hiroki Funashima in Kobe, 30 May 2019
#   modified by Hiroki Funashima in Kobe, 3 June 2019
#

from Metis.Espresso.QECheckFileType import QECheckFileType
from Metis.Espresso.QEin2Spg import QEin2Spg
from Metis.Espresso.QEout2Spg import QEout2Spg
from Metis.SpaceGroup.SpaceGroup import SpaceGroup
import os


class QE2Spg(object):
    def __init__(self, filename):
        self.filename = filename
        filecheck = QECheckFileType(filename)
        self.filetype = filecheck.filetype
        self.main()

    def main(self):
        name, ext = os.path.splitext(os.path.basename(self.filename))
        configfile = name + '.crystal_structure.in'
        if self.filetype == 'inputfile':
            QEin2Spg(self.filename, configfile)
        elif self.filetype == 'optimize':
            QEout2Spg(self.filename, configfile)
        elif self.filename is None:
            return None
        self.space_group = SpaceGroup(configfile)

