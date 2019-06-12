#!/usr/bin/env python3
import glob
import os
import subprocess

for ciffile in glob.glob('*.cif'):
    command = '../../convert_cif2spg.py {}'.format(ciffile)
    print('test ...{}'.format(ciffile))
    proc = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)
    stdout, stderr = proc.communicate()
    name, ext = os.path.splitext(os.path.basename(ciffile))
    outfile = '{}.symmetry.out'.format(name)
    with open(outfile, 'w') as fout:
        fout.write(stdout.decode('utf-8'))
    print('-> intemediate file:{}.crystal_structure.in'.format(name))
    print('-> result:{}'.format(outfile))
    print()
