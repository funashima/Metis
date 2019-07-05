#!/usr/bin/env python3
from Metis.Base.ParseConfig import ParseConfig
from Metis.Base.TspaceToolbox import TspaceToolbox
import pymatgen as mg


class GenCif(TspaceToolbox):
    def __init__(self, crystal_data):
        self.configure = ParseConfig(crystal_data)
        self.main()

    def main(self, filename='crystal.cif'):
        a, b, c = self.configure.lattice_length[:3]
        a, b, c = [self.bohr2ang(x) for x in [a, b, c]]
        alpha, beta, gamma = self.configure.lattice_angle[:3]
        lattice = mg.Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        atoms = self.configure.atom_list
        atomic_position = self.configure.atomic_position
        structure = mg.Structure(lattice, atoms, atomic_position)
        structure.to(filename=filename)
        print('generated:{}'.format(filename))
