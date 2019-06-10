#!/usr/bin/env python3
import os


class ParseGenerator(object):
    def __init__(self, generator_file):
        if os.path.isfile(generator_file):
            self.generator_file = generator_file
        else:
            print('file:{} is not found.'.format(generator_file))
            exit()
        self.parse()
        self.data_check()

    def parse(self):
        self.generator_list = [[] for x in range(230)]
        for line in open(self.generator_file):
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
            self.ispg, self.il, self.ngen = list(map(lambda x: int(x),
                                                     data[:3]))
            self.schname, self.hmname = data[3:5]
        elif data_type == 'operator':
            rot = int(data[0])
            tvec = list(map(lambda x: int(x), data[1:7]))
            self.generator.append({'rotation': rot, 'translation': tvec})
        else:
            print('unknown data_type')
            exit()


if __name__ == '__main__':
    obj = ParseGenerator('generator')
    obj.show_info()
