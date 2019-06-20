#!/usr/bin/env python3
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
import sys

def get_args():
    if len(sys.argv) < 2:
        print('usage {0} [ispg] [ichoice]'.format(sys.argv[0]))
        exit()
    ispg = sys.argv[1]
    if len(sys.argv) == 2:
        ichoice = 1
    else:
        ichoice = int(sys.argv[2])
    return [ispg, ichoice]


if __name__ == '__main__':
    ispg, ichoice = get_args()
    obj = ParseGenerator('generator')
    for gen in obj.get_generator(ispg, ichoice=ichoice):
        print(gen)
