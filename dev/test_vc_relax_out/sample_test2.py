#!/usr/bin/env python3
import glob
import os
import subprocess

for ciffile in glob.glob('*crystal*.in'):
    #command = '../QEout2Spg.py {}'.format(ciffile)
    command = '../sample1.py {}'.format(ciffile)
    print('test ...{}'.format(ciffile))
    proc = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)
    stdout, stderr = proc.communicate()
    print(stdout.decode('utf-8'))
