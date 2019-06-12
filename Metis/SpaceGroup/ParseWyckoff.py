#!/usr/bin/env python3
#
import os


class ParseWyckoff(object):
    def __init__(self, wyckoff_data):
        if os.path.isfile(wyckoff_data):
            self.wyckoff_data = wyckoff_data
        else:
            print('file:{} is not found.'.format(wyckoff_data))
            exit()
        self.parse()
        self.data_check()

    def show_info(self):
        for i in range(230):
            print('space group = {}'.format(i+1))
            for (ichoice, pos) in enumerate(self.wyckoff_position[i]):
                print('  index of choice = {}'.format(ichoice+1))
                for site_info in pos.position:
                    print('    site name: {}'.format(site_info['site_name']))

    def parse(self):
        self.wyckoff_position = [[] for x in range(230)]
        ispg = None
        for line in open(self.wyckoff_data, 'r'):
            linebuf = line.strip()
            for comment_symbol in ['#', '!']:
                if comment_symbol in linebuf:
                    linebuf = linebuf.split(comment_symbol)[0].strip()
            if linebuf == '':
                continue
            if len(linebuf.split()) == 2:
                ispg = int(linebuf.split()[0])
                wyckoff_obj = self.wyckoff_position[ispg-1]
                wyckoff_obj.append(WyckoffPositionData())
                wyckoff_obj[-1].set_data(linebuf, data_type='number')
            elif len(linebuf.split()) == 4:
                if ispg is None:
                    print('Index of Space Group is not defined')
                    exit()
                wyckoff_obj = self.wyckoff_position[ispg-1]
                wyckoff_obj[-1].set_data(linebuf, data_type='position')
            else:
                print('data format is incorrect.')
                exit()

    def data_check(self):
        for ispg in range(230):
            for (ichoice, pos) in enumerate(self.wyckoff_position[ispg]):
                if not pos.npos == len(pos.position):
                    print('number of atmic position is incompatible')
                    print('space group = {}'.format(ispg + 1))
                    print('index of choice = {}'.format(ichoice + 1))
                    print('npos = {}'.format(pos.npos))
                    print('num of op = {}'.format(len(pos.position)))
                    exit()


class WyckoffPositionData(object):
    def __init__(self):
        self.position = []

    def set_data(self, linebuf, data_type):
        data = linebuf.split()
        if data_type == 'number':
            self.ispg, self.npos = [int(x) for x in data[:2]]
        elif data_type == 'position':
            if not len(data) == 4:
                print('data format is incorrect in wyckoff file(position)')
                print('> {}'.format(linebuf))
                exit()
            site_name = data[0]
            wyckoff_position = data[1:4]
            position = {'site_name': site_name, 'position': wyckoff_position}
            self.position.append(position)
        else:
            print('data format is incorrect in wyckoff file')
            print('> {0} : space group = {1}'.format(linebuf, self.ispg))
            exit()


if __name__ == '__main__':
    obj = ParseWyckoff('wycoff')
    obj.show_info()
