#!/usr/bin/env python
import glob
import os

for target in ['*crystal_structure.in', '*.symmetry.out']:
    for filename in glob.glob(target):
        print('remove -> {}'.format(filename))
        os.remove(filename)
