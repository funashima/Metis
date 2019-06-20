#!/usr/bin/env python3
import os
home = os.environ['HOME']
pslib = os.path.join(home, 'QE/pseudo/pslibrary.1.0.0/')
pslib_version = '1.0.0'

class SelectPseudoPotential(object):
    def __init__(self, pslib, pslib_version):
        if os.path.isdir(pslib):
            self.pslib = pslib
            self.pslib_version = pslib_version
        else:
            print('dir:{} is not found.'.format(pslib))
            exit()
        
    def get_pseudo_potential(self, element, dft_type='pbesol',
                             semi_core=True, pp_type='USPP'):
        if semi_core:
            sc_type = 'dn'
        else:
            sc_type = 'n'
        if dft_type == 'pbesol':
            pre_dirname = dft_type
        elif dft_type == 'pz':
            pre_dirname = dft_type
        elif dft_type == 'pw91':
            pre_dirname = dft_type
        elif dft_type == 'bp':
            pre_dirname = dft_type
        elif dft_type == 'rel-bp':
            pre_dirname = dft_type
        elif dft_type == 'rel-pbe':
            pre_dirname = dft_type

        dirname = self.pslib + '{}/PSEUDOPOTENTIALS/'.format(pre_dirname)

        if pp_type == 'USPP':
            pp_name = 'rrkjus'
        elif pp_type == 'PAW':
            pp_name = 'kjpaw'
        else:
            print('unknown pp_type')
            exit()

        ps_name = '{0}.{1}-{2}-{3}_psl.{4}.UPF'.format(element, dft_type,
                                                      sc_type, pp_name,
                                                      self.pslib_version)
        pp_file = '{0}{1}'.format(dirname, ps_name)
        if os.path.isfile(pp_file):
            return pp_file
        else:
            return None

if __name__ == '__main__':
    element = 'Se'
    semi_core = True
    #dft_type = 'pbesol'
    #dft_type = 'pz'
    dft_type = 'pw91'
    pp_type = 'PAW'
    pp = SelectPseudoPotential(pslib, pslib_version)

    for element in ['Se', 'Te', 'Bi']:
        ppfile = pp.get_pseudo_potential(element,
                                         semi_core=True,
                                         pp_type=pp_type)
        print(ppfile)

        

