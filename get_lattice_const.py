#!/usr/bin/env python3
#
# ref https://docs.python.org/ja/3/library/random.html
#
# generate lattice constant, a, b, c, alpha, beta and gamma automatically
# written by Hiroki Funashima in Kobe, 4 June 2019
#

from Metis.Structure.GetRandomLatticeConstant import GetRandomLatticeConstant


if __name__ == '__main__':
    const_volume = 1000
    max_coa_ratio = 4.0
    #  crystal_system = 'ortho'
    crystal_system = 'triclinic'
    #  crystal_system = 'hexagonal'
    lattice = GetRandomLatticeConstant(const_volume=const_volume,
                                       crystal_system=crystal_system)

    #  variate = 'random'
    #  variate = 'normal'
    max_trial = 10
    #  print('variate:{}'.format(variate))
    for x in range(max_trial):
        print('*** trial random length index = {}'.format(x))
        a, b, c = lattice.get_lattice_length(max_coa_ratio=max_coa_ratio)
        print('  a = {:>6.2f}'.format(a))
        print('  b = {:>6.2f}'.format(b))
        print('  c = {:>6.2f}'.format(c))
        print('c/a = {:>6.2f}'.format(c/a))
        print('--')
        alpha, beta, gamma = lattice.get_lattice_angle()
        print('alpha = {:>6.2f}'.format(alpha))
        print('beta  = {:>6.2f}'.format(beta))
        print('gamma = {:>6.2f}'.format(gamma))
        print('--')
        new_vol = lattice.volume(a, b, c)
        err = abs(const_volume - new_vol)
        print('  vol = {0:5.2f} (err={1:5.2e})'.format(new_vol, err))
        err = abs(const_volume - lattice.volume(a, b, c))
        print()
    print()
