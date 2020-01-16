#!/usr/bin/env python3
import os
import glob

def check_converged(filename='espresso_relax.out'):
    if not os.path.isfile(filename):
        print('file: {} is not found'.format(filename))
        exit()
    for line in open(filename, 'r'):
        linebuf = line.strip()
        if 'final enthalpy' in linebuf.lower():
            return True
    return False

for dirname in glob.glob('Si*'):
    filename = os.path.join(dirname, 'espresso_relax.out')
    if os.path.isfile(filename):
        print('{} ... '.format(dirname), end='')
        if check_converged(filename):
            print('ok')
        else:
            print('not converged')

#print(check_converged())


