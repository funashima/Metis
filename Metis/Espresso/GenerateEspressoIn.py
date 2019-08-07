#!/usr/bin/env python3
#
#

import os
import re
import shutil
import pymatgen as mg
from Metis.Structure.GenerateCrystal import GenerateCrystal
from Metis.Structure.ParseConfigStructure import ParseConfigStructure
from Metis.Espresso.SelectPseudoPotential import SelectPseudoPotential
from Metis.Espresso.ParseQEPseudoPotential import ParseQEPseudoPotential
from Metis.Espresso.ParseConfigQE import ParseConfigQE
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.Espresso.GenerateJobScript import GenerateJobScript
from Metis.Espresso.CheckRedundancy import CheckRedundancy
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
import subprocess


class GenerateEspressoIn(object):
    def __init__(self, configfile,
                 qe_inputfile='espresso_relax.in',
                 ispg=None,
                 sub_index=None,
                 atom_info=None,
                 submit_job=False,
                 logfile=None):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('===== Error(GenerateEspressoIn) ======')
            print('file:{} is not found.'.format(configfile))
            exit()
        #
        # temporary dir:
        #  if you need, you can change it.
        #
        self.tmpdir = os.path.join('/work', os.environ['USER'])

        self.ispg = ispg
        self.sub_index = sub_index
        self.atom_info = atom_info
        self.qe_inputfile = qe_inputfile
        self.logfile = logfile
        self.main()

        if self.eliminate_redundancy:
            compound_name = self.crystal_structure.compound_name
            redundancy = CheckRedundancy()
            consistency, prim_cell = redundancy.check(compound_name,
                                                      self.ispg,
                                                      self.sub_index)
            spg_obj = ParseGenerator()
            hmname_try = spg_obj.get_hmname(self.ispg)
            if not consistency:
                if prim_cell.ispg is None:
                    hmname_true = hmname_try
                else:
                    hmname_true = spg_obj.get_hmname(prim_cell.ispg)
                if hmname_try == hmname_true:
                    fout = open(self.logfile, 'a')
                    if compound_name == prim_cell.compound_name:
                        consistency = True
                        print(' Acceptation: redundancy check was passed ',
                              end='')
                        print('for DIR:{}.'.
                              format(self.wkdir), end='')
                        print(' This calculation will be performed.')
                        fout.write(' Acceptation: ')
                        fout.write('redundancy check was passed ')
                        fout.write('for DIR:{}.'.
                                   format(self.wkdir))
                        fout.write(' This calculation will be performed.\n')
                        fout.close()
                    else:
                        print(' Redundancy : For {0:6s},'.
                              format(self.wkdir), end='')
                        fout.write(' Redundancy : For {0:6s},'.
                                   format(self.wkdir))
                        print(' cell is over size', end='')
                        fout.write(' cell is over size')
                        print(' Alternative calculation ', end='')
                        print('will be performed in {}.'.
                              format(prim_cell.dirname))
                        fout.write(' Alternative calculation ')
                        fout.write('will be performed in {}.\n'.
                                   format(prim_cell.dirname))
                        fout.close()
                        for (i, atom) in enumerate(prim_cell.atom_info):
                            if atom['element'] != self.atom_info[i]['element']:
                                print('=== Error(GenerateEspressoIn) ===')
                                print('logic error for prim_cell.atom_info')
                                exit()
                            prim_cell.atom_info[i]['semi_core'] =\
                                self.atom_info[i]['semi_core']
                        GenerateEspressoIn(configfile=self.configfile,
                                           ispg=prim_cell.ispg,
                                           qe_inputfile=qe_inputfile,
                                           sub_index=prim_cell.sub_index,
                                           atom_info=prim_cell.atom_info,
                                           submit_job=submit_job,
                                           logfile=self.logfile)
                else:  # space group is different case
                    fout = open(self.logfile, 'a')
                    print(' Redundancy : For {0:6s},'.
                          format(self.wkdir), end='')
                    fout.write(' Redundancy : For {0:6s},'.
                               format(self.wkdir))
                    print(' you assumed the symmetry:{0}, '.
                          format(hmname_try), end='')
                    print('but metis found higher symmetry:{0}.'.
                          format(hmname_true), end='')
                    fout.write(' you assumed the symmetry:{0}, '.
                               format(hmname_try))
                    fout.write('but metis found higher symmetry:{0}.'.
                               format(hmname_true))
                    for (i, atom) in enumerate(prim_cell.atom_info):
                        if atom['element'] != self.atom_info[i]['element']:
                            print('=== Error(GenerateEspressoIn) ===')
                            print('logic error for prim_cell.atom_info')
                            exit()
                        prim_cell.atom_info[i]['semi_core'] =\
                            self.atom_info[i]['semi_core']

                    if compound_name == prim_cell.compound_name:
                        print(' This calculation will be skipped.')
                        fout.write(' This calculation will be skipped.\n')
                        fout.close()
                    else:
                        if os.path.isdir(prim_cell.dirname):
                            print(' This calculation will be skipped.')
                            fout.write(' This calculation will be skipped.\n')
                            fout.close()
                        else:
                            print(' Alternative calculation ', end='')
                            print('will be performed in {}.'.
                                  format(prim_cell.dirname))
                            fout.write(' Alternative calculation ')
                            fout.write('will be performed in {}.\n'.
                                       format(prim_cell.dirname))
                            fout.close()
                            os.makedirs(prim_cell.dirname)
                            qe_inputfile = os.path.basename(self.qe_inputfile)
                            GenerateEspressoIn(configfile=self.configfile,
                                               ispg=prim_cell.ispg,
                                               qe_inputfile=qe_inputfile,
                                               sub_index=prim_cell.sub_index,
                                               atom_info=prim_cell.atom_info,
                                               submit_job=submit_job,
                                               logfile=self.logfile)
                if not consistency:
                    if os.path.isdir(self.wkdir):
                        shutil.rmtree(self.wkdir)
                    return
            else:
                fout = open(self.logfile, 'a')
                print(' Acceptation: redundancy check was passed ', end='')
                print('for DIR:{}.'.
                      format(self.wkdir), end='')
                print(' This calculation will be performed.')
                fout.write(' Acceptation: redundancy check was passed ')
                fout.write('for DIR:{}.'.
                           format(self.wkdir))
                fout.write(' This calculation will be performed.\n')
                fout.close()

        if submit_job:
            self.submit_job()

    def main(self):
        self.set_crystal_structure()
        self.get_pseudo_potential_files()
        self.generate_inputfile()
        self.copy_pp_files()
        self.gen_script()

    def set_crystal_structure(self):
        input_data = ParseConfigStructure(self.configfile)
        self.crystal_structure = GenerateCrystal(ispg=self.ispg,
                                                 ichoice=1,
                                                 max_coa_ratio=input_data.
                                                 max_coa_ratio,
                                                 apf=input_data.apf,
                                                 delta_apf=input_data.
                                                 delta_apf,
                                                 max_try=input_data.
                                                 max_try,
                                                 thr_bond_ratio=input_data.
                                                 thr_bond_ratio,
                                                 atom_info=self.atom_info,
                                                 progress=False)
        self.set_working_area()
        self.eliminate_redundancy = input_data.eliminate_redundancy

    def set_working_area(self):
        compound_name = self.crystal_structure.compound_name
        prefix = '{0}_{1}_{2}'.format(compound_name, self.ispg, self.sub_index)

        #
        # set calculation directory
        #
        self.wkdir = prefix
        if os.path.isdir(self.wkdir):
            shutil.rmtree(self.wkdir)
        os.makedirs(self.wkdir)

        #
        # set pseudo potential dir
        #
        self.pseudo_dir = os.path.join(self.wkdir, 'pseudo')
        if os.path.isdir(self.pseudo_dir):
            shutil.rmtree(self.pseudo_dir)
        os.makedirs(self.pseudo_dir)

        #
        # set output dir
        #
        output_dir = os.path.join(self.wkdir, 'out')
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)

        #
        # set temporary dir
        #
        self.wfcdir = os.path.join(self.tmpdir,
                                   os.path.basename(self.wkdir))
        if os.path.isdir(self.wfcdir):
            shutil.rmtree(self.wfcdir)
        os.makedirs(self.wfcdir)

        #
        # set name about quantum espresso inputfile
        #
        self.qe_inputfile = os.path.join(self.wkdir, self.qe_inputfile)

        #
        # generate data file for crystal structure
        #
        self.crystal_structure_out = os.path.join(self.wkdir,
                                                  'crystal_structure.out')
        self.crystal_structure.show_info(filename=self.crystal_structure_out)

        #
        # infomation for psuedo potential data
        #
        self.pseudo_potential_out = os.path.join(self.wkdir,
                                                 'pseudo_potential.out')
        if os.path.isfile(self.pseudo_potential_out):
            os.remove(self.pseudo_potential_out)

    def copy_pp_files(self):
        for source_pp_file in self.pseudo_potential_file_list:
            ppot = os.path.basename(source_pp_file)
            my_pp_file = os.path.join(self.pseudo_dir, ppot)
            if os.path.isfile(source_pp_file):
                if os.path.isfile(my_pp_file):
                    os.remove(my_pp_file)
                shutil.copy2(source_pp_file, my_pp_file)

    def gen_script(self):
        GenerateJobScript(configfile=self.configfile,
                          wkdir=self.wfcdir,
                          calc_dir=self.wkdir,
                          inputfile=os.path.basename(self.qe_inputfile))

    def submit_job(self):
        command = 'cd {0} ; qsub job.py'.format(self.wkdir)
        subprocess.call(command, shell=True)

    def get_pseudo_potential_files(self):
        configure = ParseConfigQE(self.configfile)
        self.is_spin_orbit = configure.spin_orbit
        configure.atom_info = self.crystal_structure.atom_info
        self.ecutwfc = configure.ecutwfc
        self.ecutrho = configure.ecutrho
        self.use_auto_cutoff = configure.use_auto_cutoff
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

        self.pseudo_potential_file_list = []
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
                self.pseudo_potential_file_list.append(pp_file)
                #
                # generate object for pseudopotential data
                #
                pseudo_potential = ParseQEPseudoPotential(pp_file)
                pseudo_potential.show_info(self.pseudo_potential_out)

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
        configure = ParseConfigQE(self.configfile)
        degauss = configure.degauss
        electron_maxstep = configure.electron_maxstep
        conv_thr = configure.conv_thr
        mixing_beta = configure.mixing_beta
        kpoints = configure.kpoints
        press = configure.press
        cell_factor = configure.cell_factor
        etot_conv_thr = configure.etot_conv_thr
        forc_conv_thr = configure.forc_conv_thr

        lattice_length = self.crystal_structure.lattice_length
        if self.crystal_structure.il == -1:
            is_rhombohedral = True
        else:
            is_rhombohedral = False
        #
        # check cutoff energy
        #
        if self.ecutwfc is None:
            if self.use_auto_cutoff:
                self.ecutwfc = self.rcutoff[0][1]
            else:
                print('===== Error(ecutwfc) =====')
                print(' ecutwfc is not defined.')
                print(' Suggested minimum cutoff ', end='')
                print('for wavefunctions: {:8.5f} (Ry)'.
                      format(self.rcutoff[0][1]))
                exit()
        if self.ecutrho is None:
            if self.use_auto_cutoff:
                self.ecutrho = self.rcutoff[1][1]
            else:
                self.ecutrho = self.ecutwfc * 12.0

        QEInWriteControl(filename=self.qe_inputfile,
                         wfcdir=self.wfcdir,
                         etot_conv_thr=etot_conv_thr,
                         forc_conv_thr=forc_conv_thr)
        QEInWriteSystem(filename=self.qe_inputfile,
                        ntypes=self.nelements,
                        natoms=self.total_atoms,
                        nbands=self.total_z,
                        ecutwfc=self.ecutwfc,
                        ecutrho=self.ecutrho,
                        degauss=degauss,
                        crystal_structure=self.crystal_structure,
                        is_spin_orbit=self.is_spin_orbit)
        QEInWriteElectrons(filename=self.qe_inputfile,
                           electron_maxstep=electron_maxstep,
                           mixing_beta=mixing_beta,
                           conv_thr=conv_thr)
        QEInWriteIons(filename=self.qe_inputfile)
        QEInWriteCell(filename=self.qe_inputfile,
                      press=press,
                      cell_factor=cell_factor)
        QEInWriteOthers(filename=self.qe_inputfile,
                        pseudo_potentials=self.ppots,
                        atom_info=self.crystal_structure.atom_info,
                        kpoints=kpoints,
                        lattice_length=lattice_length,
                        is_rhombohedral=is_rhombohedral)


class QEInWriteControl(object):
    def __init__(self, filename='espresso_relax.in',
                 wfcdir=None, etot_conv_thr=1.0e-4, forc_conv_thr=1.0e-3):
        self.filename = filename
        self.wfcdir = wfcdir
        self.etot_conv_thr = etot_conv_thr
        self.forc_conv_thr = forc_conv_thr
        self.main()

    def main(self):
        indent = '   '
        with open(self.filename, 'w') as fout:
            fout.write(' &control\n')
            fout.write("{0}calculation = 'vc-relax',\n".format(indent))
            fout.write("{0}restart_mode = 'from_scratch',\n".format(indent))
            fout.write("{0}pseudo_dir = './pseudo/',\n".format(indent))
            fout.write("{0}outdir = './out/',\n".format(indent))
            fout.write("{0}wfcdir = '{1}',\n".format(indent, self.wfcdir))
            fout.write("{0}disk_io = 'minimal',\n".format(indent))
            fout.write("{0}etot_conv_thr = {1:5.3e},\n".
                       format(indent, self.etot_conv_thr))
            fout.write("{0}forc_conv_thr = {1:5.3e},\n".
                       format(indent, self.forc_conv_thr))
            fout.write("{0}tstress = .true.,\n".format(indent))
            fout.write("{0}tprnfor = .true.,\n".format(indent))
            fout.write(' /\n\n')


class QEInWriteSystem(TspaceToolbox):
    def __init__(self, filename='espresso_relax.in',
                 ntypes=None, natoms=None, nbands=None,
                 ecutwfc=None, ecutrho=None,
                 degauss=0.2,
                 crystal_structure=None,
                 is_spin_orbit=False):
        self.filename = filename
        self.ntypes = ntypes
        self.natoms = natoms
        self.nbands = nbands
        self.ecutwfc = ecutwfc
        self.ecutrho = ecutrho
        self.degauss = degauss
        self.crystal_structure = crystal_structure
        self.is_spin_orbit = is_spin_orbit
        self.get_ibrav()
        self.main()

    def get_ibrav(self):
        self.celldm = [None] * 6
        self.crystal_system = self.crystal_structure.crystal_system
        self.il = self.crystal_structure.il
        if re.search('^cub', self.crystal_system):
            self.set_cubic()
        elif re.search('^tet', self.crystal_system):
            self.set_tetragonal()
        elif re.search('^ortho', self.crystal_system):
            self.set_orthorhombic()
        elif re.search('^monoc', self.crystal_system):
            self.set_monoclinic()
        elif re.search('^tric', self.crystal_system):
            self.set_triclinic()
        elif re.search('^hex', self.crystal_system):
            self.set_hexagonal()
        elif re.search('^rhombo', self.crystal_system):
            self.set_rhombohedral()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &system\n')
            fout.write('   ibrav = {},\n'.format(self.ibrav))
            fout.write('   celldm(1) = {:10.5f},\n'.format(self.celldm[0]))
            for i in range(1, 6):
                if self.celldm[i] is not None:
                    fout.write('   celldm({0}) = {1:10.5f},\n'
                               .format(i+1, self.celldm[i]))
            fout.write('   nat = {0},'.format(self.natoms))
            fout.write(' ntyp = {0},'.format(self.ntypes))
            fout.write(' nbnd = {0},\n'.format(round(self.nbands)))
            fout.write("   occupations='smearing', degauss={:.3g},"
                       .format(self.degauss))
            fout.write(" smearing='mp',\n")
            fout.write('   ecutwfc = {0:8.4f},'.format(self.ecutwfc))
            fout.write('   ecutrho = {0:8.4f},\n'.format(self.ecutrho))
            if self.is_spin_orbit:
                fout.write('   lspinorb = .true.,\n')
                fout.write('   noncolin = .true.,\n')
            fout.write(' /\n\n')

    def set_cubic(self):
        a = self.crystal_structure.lattice_length[0]
        self.celldm[0] = self.ang2bohr(a)
        if self.il == 1:
            self.ibrav = 1
        elif self.il == 2:
            self.ibrav = 2
        elif self.il == 3:
            self.ibrav = 3

    def set_tetragonal(self):
        a = self.crystal_structure.lattice_length[0]
        c = self.crystal_structure.lattice_length[2]
        self.celldm[0] = self.ang2bohr(a)
        self.celldm[2] = c/a
        if self.il == 1:
            self.ibrav = 6
        if self.il == 3:
            self.ibrav = 7

    def set_orthorhombic(self):
        a = self.crystal_structure.lattice_length[0]
        b = self.crystal_structure.lattice_length[1]
        c = self.crystal_structure.lattice_length[2]
        self.celldm[0] = self.ang2bohr(a)
        self.celldm[1] = b/a
        self.celldm[2] = c/a
        if self.il == 1:
            self.ibrav = 8
        elif self.il == 4:
            self.ibrav = 9
        elif self.il == 2:
            self.ibrav = 10
        elif self.il == 3:
            self.ibrav = 11

    def set_monoclinic(self):
        a = self.crystal_structure.lattice_length[0]
        b = self.crystal_structure.lattice_length[1]
        c = self.crystal_structure.lattice_length[2]
        gamma = self.crystal_structure.lattice_angle[2]
        self.celldm[0] = self.ang2bohr(a)
        self.celldm[1] = b/a
        self.celldm[2] = c/a
        self.celldm[3] = self.deg2cosine(gamma)
        if self.il == 1:
            self.ibrav = 12
        elif self.il == 4:
            self.ibrav = 13

    def set_triclinic(self):
        a = self.crystal_structure.lattice_length[0]
        b = self.crystal_structure.lattice_length[1]
        c = self.crystal_structure.lattice_length[2]
        alpha = self.crystal_stucture.lattice_angle[0]
        beta = self.crystal_stucture.lattice_angle[1]
        gamma = self.crystal_stucture.lattice_angle[2]
        self.celldm[0] = self.ang2bohr(a)
        self.celldm[1] = b/a
        self.celldm[2] = c/a
        self.celldm[3] = self.deg2cosine(alpha)
        self.celldm[4] = self.deg2cosine(beta)
        self.celldm[5] = self.deg2cosine(gamma)
        self.ibrav = 14

    def set_hexagonal(self):
        a = self.crystal_structure.lattice_length[0]
        c = self.crystal_structure.lattice_length[2]
        self.celldm[0] = self.ang2bohr(a)
        self.celldm[2] = c/a
        self.ibrav = 4

    def set_rhombohedral(self):
        a = self.crystal_structure.lattice_length[0]
        c = self.crystal_structure.lattice_length[2]
        a_trg, alpha_trg = self.hex2trig_lattice_params(a, c)
        self.celldm[0] = self.ang2bohr(a_trg)
        self.celldm[3] = self.deg2cosine(alpha_trg)
        self.ibrav = 5


class QEInWriteElectrons(object):
    def __init__(self, filename='espresso_relax.in',
                 electron_maxstep=100,
                 mixing_beta=0.7,
                 conv_thr=1.0e-6):
        self.electron_maxstep = electron_maxstep
        self.mixing_beta = mixing_beta
        self.conv_thr = conv_thr
        self.filename = filename
        self.main()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &electrons\n')
            fout.write("   diagonalization = 'david',\n")
            fout.write("   electron_maxstep = {},\n"
                       .format(self.electron_maxstep))
            fout.write("   mixing_mode = 'plain',\n")
            fout.write("   mixing_beta = {:.5g},\n"
                       .format(self.mixing_beta))
            fout.write("   conv_thr = {:8.5e},\n"
                       .format(self.conv_thr))
            fout.write(' /\n\n')


class QEInWriteIons(object):
    def __init__(self, filename='espresso_relax.in'):
        self.filename = filename
        self.main()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &ions\n')
            fout.write("   ion_dynamics = 'bfgs',\n")
            fout.write(' /\n\n')


class QEInWriteCell(object):
    def __init__(self, filename='espresso_relax.in',
                 press=0.0, cell_factor=1.0):
        self.filename = filename
        self.press = press
        self.cell_factor = cell_factor
        self.main()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write(' &cell\n')
            fout.write("   cell_dynamics = 'bfgs',\n")
            fout.write("   press = {:.3f}, \n".format(self.press))
            fout.write("   cell_factor = {:.3f},\n".format(self.cell_factor))
            fout.write(' /\n\n')


class QEInWriteOthers(TspaceToolbox):
    def __init__(self, filename='espresso_relax.in',
                 pseudo_potentials=None, atom_info=None,
                 kpoints=[6, 6, 6], lattice_length=None,
                 is_rhombohedral=False):
        self.filename = filename
        self.ppots = pseudo_potentials
        self.atom_info = atom_info
        self.kpoints = kpoints
        self.lattice_length = lattice_length
        self.is_rhombohedral = is_rhombohedral
        self.main()

    def main(self):
        with open(self.filename, 'a') as fout:
            fout.write('ATOMIC_SPECIES\n')
            for ppot in self.ppots:
                element_symbol = ppot.split('.')[0]
                mass = mg.Element(element_symbol).atomic_mass
                fout.write(' {0:>2s} '.format(element_symbol))
                fout.write('  {0:8.4f}'.format(mass))
                fout.write('  {}\n'.format(ppot))
            fout.write('\n')

            fout.write('ATOMIC_POSITIONS {crystal}\n')
            for atom in self.atom_info:
                element = atom['element']
                for i in range(len(atom['wyckoff_position'])):
                    for r in atom['positions'][i]:
                        if self.is_rhombohedral:
                            pos = self.hex2trig(r)
                        else:
                            pos = r
                        fout.write(' {:>2s}'.format(element))
                        for j in range(3):
                            fout.write('  {:9.6f}'.
                                       format(float(pos[j])))
                        fout.write('\n')
            fout.write('\n')
            fout.write('K_POINTS {automatic}\n')
            for i in range(3):
                fout.write(' {}'.format(self.kpoints[i]))
            fout.write(' 0 0 0\n')
