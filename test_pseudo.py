#!/usr/bin/env python3
import os
import glob


class SelectPseudoPotential(object):
    def __init__(self, pslib):
        if os.path.isdir(pslib):
            self.pslib = pslib
        else:
            print('dir:{} is not found.'.format(pslib))
            exit()
        
    def get_pseudo_potential(self, element,
                             dft_type='pbesol',
                             spin_orbit=False,
                             semi_core=[], pp_type='USPP'):
        electron_config = 'n'
        for x in semi_core:
            for lname in ['s', 'p', 'd', 'f']:
                if x == lname:
                    electron_config = lname + electron_config

        if dft_type == 'pbe':
            pre_dirname = dft_type
        elif dft_type == 'pbesol':
            pre_dirname = dft_type
        elif dft_type == 'pz':
            pre_dirname = dft_type
        elif dft_type == 'pw91':
            pre_dirname = dft_type
        elif dft_type == 'bp':
            pre_dirname = dft_type

        if spin_orbit:
            dft_type = 'rel-' + dft_type
            pre_dirname = 'rel-' + pre_dirname

        dirname = os.path.join(self.pslib, pre_dirname, 'PSEUDOPOTENTIALS')

        if pp_type == 'USPP':
            pp_name = 'rrkjus'
        elif pp_type == 'PAW':
            pp_name = 'kjpaw'
        else:
            print('unknown pp_type')
            exit()

        ps_name = '{0}.{1}-{2}-{3}_psl.*.UPF'.format(element, dft_type,
                                                      electron_config, pp_name)
        pp_file = '{0}/{1}'.format(dirname, ps_name)
        if len(glob.glob(pp_file)) > 0:
            return glob.glob(pp_file)[0]
        return None

if __name__ == '__main__':
    pslib = '/home/funashima/QE/pseudo/pslibrary'
    dft_type = 'pbe'
    #dft_type = 'pbesol'
    #dft_type = 'pz'

    pp_type = 'PAW'
    #pp_type = 'USPP'

    spin_orbit = True

    atoms = [{'element':'Se', 'semi_core': []},
             {'element':'Te', 'semi_core': ['d']},
             {'element':'Bi', 'semi_core': ['d']}]

    pp = SelectPseudoPotential(pslib)

    for atom_info in atoms:
        element = atom_info['element']
        semi_core = atom_info['semi_core']

        ppfile = pp.get_pseudo_potential(element,
                                         semi_core=semi_core,
                                         spin_orbit=spin_orbit,
                                         pp_type=pp_type)
        print(ppfile)
