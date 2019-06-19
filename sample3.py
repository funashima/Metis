#!/usr/bin/env python3

from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff


if __name__ == '__main__':
    ispg = 21
    ichoice = 1
    wletter = 'c'
    wyckoff = ParseWyckoff()
    for wlet in ['c', 'd', 'e', 'f', 'g', 'h']:
        pos = wyckoff.get_atomic_position(ispg=ispg,
                                          ichoice=ichoice,
                                          x=None, y=None, z=None,
                                          wyckoff_letter=wlet)
        print(pos)
