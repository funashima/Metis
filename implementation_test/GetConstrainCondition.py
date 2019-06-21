#!/usr/bin/env python3
#
# get constraint condtion from QEstyle, crystal_in and so on.
#  written by Hiroki Funashima 12 June 2019 in Kobe
#
import math
from Metis.Base.ParseConfig import ParseConfig
import re
import pymatgen as mg


class GetConstrainConditionCrystalIn(object):
    def __init__(self, inputfile, unit='angstrom'):
        self.config_obj = ParseConfig(inputfile)
        self.get_volume(unit)

    def sphere_volume(self, r):
        return (4.0 / 3.0) * math.pi * (r**3)

    def get_volume(self, unit='angstrom'):
        self.volume = 0.0
        for atom_name in self.config_obj.atom_list:
            atomic_radius = mg.Element(atom_name).atomic_radius
            self.volume += self.sphere_volume(atomic_radius)
        #
        # -Atomic Packing Factor(APF)-
        #
        # APF: bcc = pi*sqrt(3.0) / 8 = 0.68
        #      fcc/hcp = pi/sqrt(18) = 0.74
        #      diamond = 0.34
        #      simple cubic lattice = 0.52
        #
        apf = 0.68
        self.volume /= apf
        if re.search('^bohr', unit):
            self.volume *= (1.8897259885789233)**3

    def ang2bohr(self, r):
        return r / 0.529177249

    def ang2aunit(self, r):
        a = self.config_obj.lattice_length[0]
        return self.ang2bohr(r) / a

    def check_bond_length(self):
        natom = self.config_obj.natom
        check = True
        for i in range(natom):
            atom1 = self.config_obj.atom_list[i]
            for j in range(i+1, natom):
                atom2 = self.config_obj.atom_list[j]
                sum_atomic_radius = mg.Element(atom1).atomic_radius +\
                    mg.Element(atom2).atomic_radius
                bond_length = self.config_obj.calc_length(i, j, aunit=True)
                if bond_length < sum_atomic_radius:
                    print('{0}-{1}:'.format(atom1, atom2), end='')
                    print('ideal = {0:8.5f} but real = {1:8.5f}'.
                          format(sum_atomic_radius, bond_length))
                    check = False
        return check
