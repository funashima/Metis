#!/usr/bin/env python3
#
# test program for GenerateCrystal
#

from Metis.Structure.GenerateCrystal import GenerateCrystal
from Metis.Structure.ParseTestInput import ParseTestInput

inputfile = 'crystal.in'
input_data = ParseTestInput(inputfile)
GenerateCrystal(ispg=input_data.space_group,
                ichoice=input_data.ichoice,
                atom_info=input_data.atom_info).show_info()
