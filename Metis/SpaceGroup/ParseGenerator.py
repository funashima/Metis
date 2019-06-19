#!/usr/bin/env python3
from Metis.SpaceGroup.LibTspaceDataBase import LibTspaceDataBase


class ParseGenerator(object):
    def __init__(self):
        self.parse()
        self.data_check()

    def parse(self):
        self.generator_list = [[] for x in range(230)]
        generator_file = LibTspaceDataBase().generator_database().split('\n')
        for line in generator_file:
            linebuf = line.strip()
            for comment_symbol in ['#', '!']:
                if comment_symbol in linebuf:
                    linebuf = linebuf.split(comment_symbol)[0].strip()
            if linebuf == '':
                continue
            ndata = len(linebuf.split())
            if ndata == 5:
                ispg = int(linebuf.split()[0])
                self.generator_list[ispg-1].append(GeneratorData())
                self.generator_list[ispg-1][-1].set_data(linebuf,
                                                         data_type='spg')
            elif ndata == 7:
                self.generator_list[ispg-1][-1].set_data(linebuf,
                                                         data_type='operator')

    def data_check(self):
        for ispg in range(230):
            for (ichoice, gen) in enumerate(self.generator_list[ispg]):
                if not gen.ngen == len(gen.generator):
                    print('number of generator is incompatible')
                    print('space group = {}'.format(ispg + 1))
                    print('index of choice = {}'.format(ichoice + 1))
                    print('ngen = {}'.format(gen.ngen))
                    print('num of op = {}'.format(len(gen.generator)))
                    print(gen)
                    exit()

    def spgname2ispg(self, spgname):
        try:
            return int(spgname)
        except ValueError:
            for i in range(230):
                for gen in self.generator_list[i]:
                    if gen.schname == spgname or gen.hmname == spgname:
                        return i+1
            return None

    def get_generator(self, spgname, ichoice=1):
        ispg = self.spgname2ispg(spgname)
        if ispg is None:
            print('space group name:{} is incorrect.'.format(spgname))
            exit()
        if len(self.generator_list[ispg - 1]) < ichoice:
            return None
        return self.generator_list[ispg - 1][ichoice - 1].generator

    def get_il(self, spgname, ichoice=1):
        ispg = self.spgname2ispg(spgname)
        if ispg is None:
            print('space group name:{} is incorrect.'.format(spgname))
            exit()
        return self.generator_list[ispg - 1][ichoice - 1].il

    def get_hmname(self, spgname, ichoice=1):
        ispg = self.spgname2ispg(spgname)
        if ispg is None:
            print('space group name:{} is incorrect.'.format(spgname))
            exit()
        return self.generator_list[ispg - 1][ichoice - 1].hmname

    def get_schname(self, spgname, ichoice=1):
        ispg = self.spgname2ispg(spgname)
        if ispg is None:
            print('space group name:{} is incorrect.'.format(spgname))
            exit()
        return self.generator_list[ispg - 1][ichoice - 1].schname

    def get_ispg(self, spgname, ichoice=1):
        return self.spgname2ispg(spgname)

    def show_info(self):
        for i in range(230):
            print('space group = {}'.format(i+1))
            for (ichoice, gen) in enumerate(self.generator_list[i]):
                print('  schname:{}'.format(gen.schname))
                print('  hmname:{}'.format(gen.hmname))
                print('  il = {}'.format(gen.il))
                print('  index of choice = {}'.format(ichoice+1))
                for gen_info in gen.generator:
                    print('    rotation: {}'.format(gen_info['rotation']),
                          end='')
                    print('  transration:{}'.format(gen_info['translation']))


class GeneratorData(object):
    def __init__(self):
        self.generator = []

    def set_data(self, linebuf, data_type):
        data = linebuf.split()
        if data_type == 'spg':
            self.ispg, self.il, self.ngen = [int(x) for x in data[:3]]
            self.schname, self.hmname = data[3:5]
        elif data_type == 'operator':
            rot = int(data[0])
            tvec = [int(x) for x in data[1:7]]
            self.generator.append({'rotation': rot, 'translation': tvec})
        else:
            print('unknown data_type')
            exit()


if __name__ == '__main__':
    obj = ParseGenerator('generator')
    obj.show_info()
