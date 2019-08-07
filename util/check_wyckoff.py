#!/usr/bin/env python3
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
import os


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('file:{} is not found.'.format(configfile))
            exit()
        self.generator = ParseGenerator()
        self.main()

    def main(self):
        self.set_default_value()
        self.parse()
        self.check_spg()

    def get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def get_key_and_value(self, linedata, token='='):
        if token not in linedata:
            return [None, None]
        key, value = [x.strip() for x in linedata.split(token)[:2]]
        key = key.lower()
        return [key, value]

    def spg2ispg(self, spg):
        return self.generator.get_ispg(spg)

    def set_default_value(self):
        self.min_spg_index = None
        self.max_spg_index = None
        self.ichoice = 1

    def check_spg(self):
        if [self.min_spg_index, self.max_spg_index] == [None] * 2:
            print('space group is not defined')
            exit()

        if self.min_spg_index is not None:
            if self.max_spg_index is None:
                self.max_spg_index = self.min_spg_index

        if self.max_spg_index is not None:
            if self.min_spg_index is None:
                self.mim_spg_index = self.max_spg_index

        if self.min_spg_index > self.max_spg_index:
            self.min_spg_index, self.max_spg_index =\
                    self.max_spg_index, self.min_spg_index

    def parse(self):
        for line in open(self.configfile, 'r'):
            linebuf = self.get_linebuf(line)
            if linebuf == '':
                continue
            if '=' in linebuf:
                key, value = self.get_key_and_value(linebuf)
                if key is None:
                    continue
                if value == '':
                    continue
                if key == 'min_spg_index':
                    self.min_spg_index = self.spg2ispg(value)
                elif key == 'max_spg_index':
                    self.max_spg_index = self.spg2ispg(value)
                elif key == 'ichoice':
                    self.ichoice = int(value)


class CheckWyckoff(object):
    def __init__(self, configfile):
        configure = ParseConfig(configfile)
        self.wyckoff = ParseWyckoff()
        self.generator = ParseGenerator()
        self.min_spg_index = configure.min_spg_index
        self.max_spg_index = configure.max_spg_index
        self.ichoice = configure.ichoice
        self.main()

    def set_spg_list(self):
        return list(range(self.min_spg_index,
                          self.max_spg_index+1))

    def reform_pos(self, pos):
        if '(' in pos and ')' in pos and '-' in pos:
             pos = pos.replace('(', '').replace(')', '')
        if ' + ' in pos:
            pos = pos.replace(' + ', '+')
        if ' - ' in pos:
            pos = pos.replace(' - ', '-')
        return pos


    def main(self):
        for spg in self.set_spg_list():
            space_group_info = GenerateSpaceGroup(spg)
            hmname = space_group_info.hmname
            schname = space_group_info.schname
            bravais_lat = space_group_info.get_bravais_lattice_name()
            print('space group = #{} '.format(spg), end='')
            print('({0}, {1})'.format(hmname, schname))
            print(' lattice:{}'.format(bravais_lat))
            wyckoff_data = self.wyckoff.\
                wyckoff_position[spg-1][self.ichoice-1].position
            for (isite, info) in enumerate(wyckoff_data):
                site_name = info['site_name']
                print(' ({:>2d}) '.format(isite+1), end='')
                print(' atomic site:{:>4s} ('.format(site_name), end='')
                pos_list = self.wyckoff.\
                    get_atomic_position(ispg=spg,
                                        ichoice=self.ichoice,
                                        wyckoff_letter=site_name)
                for (i, pos) in enumerate(pos_list):
                    print('{:>6s}'.format(self.reform_pos(pos)), end='')
                    if i == 2:
                        print(')')
                    else:
                        print(', ', end='')


if __name__ == '__main__':
    configfile = 'space_group.in'
    CheckWyckoff(configfile)
