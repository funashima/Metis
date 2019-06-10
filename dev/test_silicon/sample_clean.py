#!/usr/bin/env python
import glob
import os

for target in ['*.out', '*crystal_structure.in']:
    for filename in glob.glob(target):
        print('remove -> {}'.format(filename))
        os.remove(filename)
