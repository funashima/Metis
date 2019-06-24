#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET


class ParseQEPseudoPotential(object):
    def __init__(self, pp_file):
        if os.path.isfile(pp_file):
            self.pp_file = pp_file
        else:
            print('file:{} is not found.'.format(pp_file))
            exit()
        self.main()

    def set_attr_str_value(self, header, keyword, vtype=None):
        exec('self.{0} = header.attrib["{0}"]'.format(keyword))

    def set_attr_int_value(self, header, keyword, vtype=None):
        exec('self.{0} = int(header.attrib["{0}"])'.format(keyword))

    def set_attr_float_value(self, header, keyword, vtype=None):
        exec('self.{0} = float(header.attrib["{0}"])'.format(keyword))

    def set_attr_bool_value(self, header, keyword, vtype=None):
        bool_value = eval('header.attrib["{0}"]'.format(keyword))
        if re.search('^(t|T)', bool_value):
            exec('self.{0} = True'.format(keyword))
        else:
            exec('self.{0} = False'.format(keyword))

    def main(self):
        tree = ET.parse(self.pp_file)
        root = tree.getroot()
        for header in root.iter('PP_HEADER'):
            for keyword in ['element', 'pseudo_type', 'relativistic']:
                self.set_attr_str_value(header, keyword)

            for keyword in ['l_max', 'l_max_rho', 'l_local',
                            'mesh_size', 'number_of_wfc', 'number_of_proj']:
                self.set_attr_int_value(header, keyword)

            for keyword in ['z_valence',
                            'wfc_cutoff', 'rho_cutoff',
                            'total_psenergy']:
                self.set_attr_float_value(header, keyword)

            for keyword in ['is_ultrasoft', 'is_paw',
                            'is_coulomb', 'has_so',
                            'has_wfc',
                            'has_gipaw', 'paw_as_gipaw',
                            'core_correction']:
                self.set_attr_bool_value(header, keyword)

    def show_info(self, filename=None):
        def sub_show_info(fout):
            fout.write('\n')
            fout.write('** Basic info:{}\n'.format(self.pp_file))
            fout.write('\n')
            fout.write('Element:{}\n'.format(self.element))
            fout.write('pseudo_type = {}\n'.format(self.pseudo_type))
            fout.write('relativistic type:{}\n'.format(self.relativistic))
            fout.write('spin-orbit ... ')
            if self.has_so:
                fout.write('yes\n')
            else:
                fout.write('no\n')
            fout.write('paw type ... ')
            if self.is_paw:
                fout.write('yes\n')
            else:
                fout.write('no\n')
            fout.write('ultrasoft type ... ')
            if self.is_ultrasoft:
                fout.write('yes\n')
            else:
                fout.write('no\n')
            fout.write('including 1/r coloumb potential ... ')
            if self.is_coulomb:
                fout.write('yes\n')
            else:
                fout.write('no\n')
            fout.write('core correction ... ')
            if self.core_correction:
                fout.write('yes\n')
            else:
                fout.write('no\n')
            fout.write('num of valence electrons = {}\n'.
                       format(self.z_valence))
            fout.write('Total Pseudo Eenrgy = {:.5g}(Ry)\n'.
                       format(self.total_psenergy))
            fout.write('cut off energy for wavefunction   = {:8.4f} (Ry)\n'.
                       format(self.wfc_cutoff))
            fout.write('cut off energy for charge density = {:8.4f} (Ry)\n'.
                       format(self.rho_cutoff))
            fout.write('\n')

        if filename is None:
            sub_show_info(sys.stdout)
        else:
            with open(filename, 'w') as fout:
                sub_show_info(fout)
