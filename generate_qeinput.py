#!/usr/bin/env python3
#
# test program for GenerateCrystal
#
import os
import re
from Metis.Structure.GenerateCrystal import GenerateCrystal
from Metis.Structure.ParseTestInput import ParseTestInput
from Metis.Espresso.SelectPseudoPotential import SelectPseudoPotential
from Metis.Espresso.ParseQEPseudoPotential import ParseQEPseudoPotential
import pymatgen as mg


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('===== Error(ParseConfig) =====')
            print('file:{} is not found.'.format(configfile))
            exit()
        self.main()

    def get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def get_key_and_value(self, linebuf):
        if '=' not in linebuf:
            return [None, None]
        key, value = [x.strip() for x in linebuf.split('=')[:2]]
        key = key.lower()
        return [key, value]

    def value2bool(self, value):
        if re.search('^(t|y)', value.lower()):
            return True
        elif re.search('^(f|n)', value.lower()):
            return False
        return False

    def check_pslib(self, value):
        if os.path.isdir(value):
            self.pslib = value
        else:
            print('===== Error(ParseConfig) =====')
            print('dir:{} is not found'.format(value))
            exit()

    def set_environments(self, linebuf):
        key, value = self.get_key_and_value(linebuf)
        if key is None:
            return
        if key == 'pslib':
            self.check_pslib(value)
        elif key == 'dft_type':
            self.dft_type = value
        elif key == 'pp_type':
            self.pp_type = value
        elif key == 'spin_orbit':
            self.spin_orbit = self.value2bool(value)

    def set_atom_info(self, linebuf):
        info = {'element': None,
                'natoms': 1,
                'semi_core': []}
        for linedata in [x.strip() for x in linebuf.split(';')]:
            key, value = self.get_key_and_value(linedata)
            if key is None:
                continue
            if key == 'element':
                info['element'] = value
            if key == 'natoms':
                info['natoms'] = int(value)
            if key == 'semi_core':
                info['semi_core'] = value.replace(',', ' ').split()
        if info['element'] is not None:
            self.atom_info.append(info)

    def set_init_value(self):
        self.dft_type = 'pbe'
        self.pp_type = 'USPP'
        self.spin_orbit = False
        self.atom_info = []

    def main(self):
        self.set_init_value()
        atom_info_region = False
        for line in open(self.configfile, 'r'):
            linebuf = self.get_linebuf(line)
            if linebuf == '':
                continue
            if atom_info_region:
                if re.search('^end_atom_info', linebuf.lower()):
                    atom_info_region = False
                    continue
                self.set_atom_info(linebuf)
            else:
                if re.search('^begin_atom_info', linebuf.lower()):
                    atom_info_region = True
                    continue
                if '=' in linebuf:
                    self.set_environments(linebuf)


class GenerateEspressoIn(object):
    def __init__(self, configfile, qe_inputfile='espresso_relax.in'):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('===== Error(GenerateEspressoIn) ======')
            print('file:{} is not found.'.format(configfile))
            exit()
        self.qe_inputfile = qe_inputfile
        self.main()

    def main(self):
        self.set_crystal_structure()
        self.get_pseudo_potential_files()
        self.generate_inputfile()

    def set_crystal_structure(self):
        input_data = ParseTestInput(self.configfile)
        self.crystal_structure = GenerateCrystal(ispg=input_data.space_group,
                                                 ichoice=input_data.ichoice,
                                                 max_coa_ratio=input_data.
                                                 max_coa_ratio,
                                                 apf=input_data.apf,
                                                 atom_info=input_data.
                                                 atom_info)

    def get_pseudo_potential_files(self):
        configure = ParseConfig(self.configfile)
        configure.atom_info = self.crystal_structure.atom_info
        pp = SelectPseudoPotential(configure.pslib)
        self.total_z = 0
        self.total_atoms = 0
        self.ppots = []
        #
        # - store list for values -
        # rcutoff = [[rmin for rwfc, rmax for rwfc],
        #            [rmin for rrho, rmax for rrho]]
        #
        self.rcutoff = [[None, None], [None, None]]

        #
        # num of kind of elements
        #
        self.nelements = len(configure.atom_info)

        for atom_info in configure.atom_info:
            element = atom_info['element']
            natoms = atom_info['natoms']
            self.total_atoms += natoms
            semi_core = atom_info['semi_core']

            pp_file = pp.get_pseudo_potential(element,
                                              dft_type=configure.dft_type,
                                              semi_core=semi_core,
                                              spin_orbit=configure.spin_orbit,
                                              pp_type=configure.pp_type)
            if pp_file is None:
                print('** Error: no pseudo potential file:{}'.format(element))
                print('   please check your condition, ', end='')
                print('specially semi_core states')
            else:
                #
                # generate object for pseudopotential data
                #
                pseudo_potential = ParseQEPseudoPotential(pp_file)

                #
                # z_valence: number of valence electron
                # rwfc: cutoff radius of wavefunction
                # rrho: cutoff radius of charge density
                #
                self.total_z += pseudo_potential.z_valence * natoms
                rwfc = pseudo_potential.wfc_cutoff
                rrho = pseudo_potential.rho_cutoff

                #
                # check minimum and maximum for rwfc and rrho
                #
                for i in range(2):
                    if self.rcutoff[0][i] is None:
                        self.rcutoff[0][i] = rwfc
                    if self.rcutoff[1][i] is None:
                        self.rcutoff[1][i] = rrho

                if self.rcutoff[0][0] > rwfc:
                    self.rcutoff[0][0] = rwfc
                if self.rcutoff[0][1] < rwfc:
                    self.rcutoff[0][1] = rwfc
                if self.rcutoff[1][0] > rrho:
                    self.rcutoff[1][0] = rrho
                if self.rcutoff[1][1] < rrho:
                    self.rcutoff[1][1] = rrho

                #
                # add pseudo potential file to list
                #
                self.ppots.append(os.path.basename(pp_file))

    def generate_inputfile(self):
        print('crystal_system = {}'.format(self.crystal_structure.crystal_system))
        print('IL = {}'.format(self.crystal_structure.il))
        print('generate:{}'.format(self.qe_inputfile))
        QEInWriteControl(filename=self.qe_inputfile)
        QEInWriteSystem(filename=self.qe_inputfile,
                        ntypes=self.nelements,
                        natoms=self.total_atoms,
                        nbands=self.total_z,
                        ecutwfc=self.rcutoff[0][1],
                        ecutrho=self.rcutoff[1][1])
        QEInWriteIons(filename=self.qe_inputfile)
        QEInWriteOthers(filename=self.qe_inputfile,
                        pseudo_potentials=self.ppots,
                        atom_info=self.crystal_structure.atom_info)


class QEInWriteControl(object):
    def __init__(self, filename='espresso_relax.in'):
        self.filename = filename
        self.main()

    def main(self):
        indent = '   '
        with open(self.filename, 'w') as fout:
            fout.write(' &control\n')
            fout.write("{0}calculation = 'vc_relax'\n".format(indent))
            fout.write("{0}restart_mode = 'from_scratch'\n".format(indent))
            fout.write("{0}pseudo_dir = './pseudo'\n".format(indent))
            fout.write("{0}diskio = 'minimal'\n".format(indent))
            fout.write("{0}tstress = .true.\n".format(indent))
            fout.write("{0}tprnfor = .true.\n".format(indent))
            fout.write(' /\n\n')


class QEInWriteSystem(object):
    def __init__(self, filename='espresso_relax.in',
                 ntypes=None, natoms=None, nbands=None,
                 ecutwfc=None, ecutrho=None, crystal_structure=None):
        self.filename = filename
        self.ntypes = ntypes
        self.natoms = natoms
        self.nbands = nbands
        self.ecutwfc = ecutwfc
        self.ecutrho = ecutrho
        self.crystal_structure = crystal_structure
        self.main()

    def get_ibrav(self, crystal_sytem=None, il=None):
        pass


    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &system\n')
            fout.write('   nat = {0},'.format(self.natoms))
            fout.write(' ntyp = {0},'.format(self.ntypes))
            fout.write(' nbnd = {0},\n'.format(round(self.nbands)))
            fout.write('   ecutwfc = {0:8.4f},'.format(self.ecutwfc))
            fout.write('   ecutrho = {0:8.4f}, \n'.format(self.ecutrho))
            fout.write(' /\n\n')


class QEInWriteIons(object):
    def __init__(self, filename='espresso_relax.in'):
        self.filename = filename
        self.main()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &ions\n')
            fout.write(' /\n\n')


class QEInWriteOthers(object):
    def __init__(self, filename='espresso_relax.in',
                 pseudo_potentials=None, atom_info=None):
        self.filename = filename
        self.ppots = pseudo_potentials
        self.atom_info = atom_info
        self.main()

    def main(self):
        #for atom in self.atom_info:
        #    print(atom)
        with open(self.filename, 'a') as fout:
            fout.write('ATOMIC_SPECIES\n')
            for ppot in self.ppots:
                element_symbol = ppot.split('.')[0]
                mass = mg.Element(element_symbol).atomic_mass
                fout.write('  {0:>2s}'.format(element_symbol))
                fout.write('  {0:8.4f}'.format(mass))
                fout.write('  {}\n'.format(ppot))
            fout.write('\n')

            fout.write('ATOMIC_POSITIONS {crystal}\n')
            for atom in self.atom_info:
                element = atom['element']
                for i in range(len(atom['wyckoff_position'])):
                    for pos in atom['positions'][i]:
                        fout.write(' {:>2s}'.format(element))
                        for j in range(3):
                            fout.write('  {:9.6f}'.
                                       format(float(pos[j])))
                        fout.write('\n')
            fout.write('\n')


if __name__ == '__main__':
    GenerateEspressoIn('crystal.in')
