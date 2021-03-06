#!/usr/bin/env python3
#
from Metis.SpaceGroup.LibTspaceDataBase import LibTspaceDataBase


class ParseWyckoff(object):
    def __init__(self):
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
        wyckoff_data = LibTspaceDataBase().wyckoff_database().split('\n')
        for line in wyckoff_data:
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

    def get_wyckoff_position(self, ispg=None, ichoice=1,
                             wyckoff_letter=None):
        if len(self.wyckoff_position[ispg - 1]) < ichoice:
            return None
        pos = self.wyckoff_position[ispg - 1][ichoice - 1].position
        input_wyckoff_letter = wyckoff_letter[-1].lower()
        for site_type in pos:
            if site_type['site_name'][-1] == input_wyckoff_letter:
                return site_type['position']
        return None

    def get_atomic_position(self, ispg=None, ichoice=1,
                            wyckoff_letter=None,
                            x=None, y=None, z=None):
        tspace_notation = self.\
                          get_wyckoff_position(ispg=ispg,
                                               ichoice=ichoice,
                                               wyckoff_letter=wyckoff_letter)
        if tspace_notation is None:
            print('====== Error(get_atomic_position in ParseWyckoff) ======')
            print(' wyckoff letter: "{}" is invalid for this space group'.
                  format(wyckoff_letter))
            print(' check international table of crystallography.')
            print()
            exit()

        tsp_table = {'0': '0',
                     'h': '1/2',
                     'q': '1/4',
                     't': '3/4',
                     '1': '1/3',
                     '2': '2/3',
                     '3': '3/8',
                     '5': '5/6',
                     '6': '1/6',
                     '7': '7/8',
                     '8': '1/8',
                     'f': '5/8',
                     'x': 'x',
                     'y': 'y',
                     'z': 'z',
                     '-': '-(x)',
                     'w': '-(y)',
                     'd': '2 * (x)',
                     'm': '-(y) + 1/2',
                     'p': 'y + 1/2',
                     'n': '-(y) + 1/4',
                     'r': 'y + 1/4',
                     's': 'x + 1/2',
                     'u': '-(x) + 1/2',
                     'v': 'x + 1/4'
                     }
        for i in range(3):
            for (keycode, op) in tsp_table.items():
                if tspace_notation[i] == keycode:
                    tspace_notation[i] = tspace_notation[i].\
                                         replace(keycode, op)
                    break
        if x is not None:
            tspace_notation = [n.replace('x', str(x)) for n in tspace_notation]
        if y is not None:
            tspace_notation = [n.replace('y', str(y)) for n in tspace_notation]
        if z is not None:
            tspace_notation = [n.replace('z', str(z)) for n in tspace_notation]
        return tspace_notation


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
    obj = ParseWyckoff()
    obj.show_info()
