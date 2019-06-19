#!/usr/bin/env python3
#
# generate LibTspaceDataBase.py
#   from tpace database files, 'wycoff' and 'generator'
# written by Hiroki Funashima in Kobe, 16 June 2019
#

database = ['./data/wycoff', './data/generator']
target = 'LibTspaceDataBase.py'
header = ['./data/header1.py', './data/header2.py']
space = ' '
space8 = ''
space12 = ''
words = ['wyckoff', 'generator']


def get_space(num):
    char = ''
    space = ' '
    for i in range(num):
        char += space
    return char


space8 = get_space(8)
space12 = get_space(12)

with open(target, 'w') as fout:
    for i in range(2):
        for line in open(header[i], 'r'):
            fout.write(line)
        for line in open(database[i], 'r'):
            fout.write('{0}{1}'.format(space12, line.rstrip() + '\n'))
        fout.write('{0}{1}\n'.format(space8, "'''"))
        fout.write('{0}return {1}\n'.format(space8, words[i]))
        if i == 0:
            fout.write('\n')
print('file:{} has been generated.'.format(target))
