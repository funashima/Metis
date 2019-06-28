#!/usr/bin/env python3
import subprocess
for i in range(3000, 4000):
    command = 'qdel {}'.format(i)
    subprocess.call(command, shell=True)
