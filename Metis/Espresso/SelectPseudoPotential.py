#!/usr/bin/env python3
#
#
#

import os
import glob


class SelectPseudoPotential(object):
    def __init__(self, pslib):
        if os.path.isdir(pslib):
            self.pslib = pslib
        else:
            print('dir:{} is not found.'.format(pslib))
            exit()

    def get_electron_config(self, semi_core):
        electron_config = 'n'
        for lname in reversed(['s', 'p', 'd', 'f']):
            for x in semi_core:
                if x == lname:
                    electron_config = lname + electron_config
        return electron_config

    def get_dirname(self, dft_type, spin_orbit):
        for exc_name in ['pbe', 'pbesol', 'pz', 'pw91', 'bp']:
            if dft_type.lower() == exc_name:
                pre_dirname = exc_name
                break
        if spin_orbit:
            pre_dirname = 'rel-' + pre_dirname
        return pre_dirname

    def get_pp_dir(self, dft_type, spin_orbit):
        pre_dirname = self.get_dirname(dft_type, spin_orbit)
        dirname = os.path.join(self.pslib,
                               pre_dirname,
                               'PSEUDOPOTENTIALS')
        return dirname

    def get_ppfile_name(self, element,
                        dft_type, spin_orbit,
                        electron_config, pp_type):
        pp_type_list = {'uspp': 'rrkjus', 'paw': 'kjpaw'}
        pp_name = pp_type_list[pp_type.lower()]
        if pp_name is None:
            print('unknown pp_type')
            exit()
        pre_dirname = self.get_dirname(dft_type, spin_orbit)
        return '{0}.{1}-{2}-{3}_psl.*.UPF'.\
               format(element, pre_dirname, electron_config, pp_name)

    def get_pseudo_potential(self, element,
                             dft_type='pbesol',
                             spin_orbit=False,
                             semi_core=[], pp_type='USPP'):
        electron_config = self.get_electron_config(semi_core)
        pp_dir = self.get_pp_dir(dft_type, spin_orbit)
        pp_filename = self.get_ppfile_name(element,
                                           dft_type, spin_orbit,
                                           electron_config, pp_type)
        abs_pp_filename = '{0}/{1}'.format(pp_dir, pp_filename)
        if len(glob.glob(abs_pp_filename)) > 0:
            return glob.glob(abs_pp_filename)[0]
        return None
