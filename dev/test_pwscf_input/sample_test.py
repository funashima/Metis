#!/usr/bin/env python3
import glob
import os
import subprocess

print('=== Test ===')
print(' test for QE2Spg.py')
print()

for ciffile in glob.glob('*.in'):
    #command = '../QE2Spg.py {}'.format(ciffile)
    command = '../../qe_analyzer.py {}'.format(ciffile)
    if 'crystal_structure' in ciffile:
        continue
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
