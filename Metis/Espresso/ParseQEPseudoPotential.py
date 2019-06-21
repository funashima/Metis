#!/usr/bin/env python3
import os
import re
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

    def show_info(self):
        print()
        print('** Basic info:{}'.format(self.pp_file))
        print()
        print('Element:{}'.format(self.element))
        print('pseudo_type = {}'.format(self.pseudo_type))
        print('relativistic type:{}'.format(self.relativistic))
        print('spin-orbit ...', end='')
        if self.has_so:
            print('yes')
        else:
            print('no')
        print('paw type ... ', end='')
        if self.is_paw:
            print('yes')
        else:
            print('no')
        print('ultrasoft type ... ', end='')
        if self.is_ultrasoft:
            print('yes')
        else:
            print('no')
        print('including 1/r coloumb potential ... ', end='')
        if self.is_coulomb:
            print('yes')
        else:
            print('no')
        print('core correction ... ', end='')
        if self.core_correction:
            print('yes')
        else:
            print('no')
        print('num of valence electrons = {}'.format(self.z_valence))
        print('Total Pseudo Eenrgy = {:.5g}(Ry)'.format(self.total_psenergy))
        print('cut off energy for wavefunction   = {:8.4f} (Ry)'.
              format(self.wfc_cutoff))
        print('cut off energy for charge density = {:8.4f} (Ry)'.
              format(self.rho_cutoff))
        print()
