#!/usr/bin/python
#
# -- job name --
#
#$ -N /work/funashima/Bi2Te3_166
#
#$ -pe fillup 1
#
#$ -V
#$ -S /usr/bin/python
#$ -cwd
#$ -e espresso_relax.err
#$ -o espresso_relax.out
#

import os
import shutil
import subprocess

wkdir = '/work/funashima/Bi2Te3_166'
bindir = '/home/funashima/QE/espresso/bin'
inputfile = 'espresso_relax.in'
if not os.path.isdir(wkdir):
    os.makedirs(wkdir, exist_ok=True)

command = "{0}/pw.x < {1}".\
    format(bindir, inputfile)
proc = subprocess.call(command, shell=True)
if proc < 1:
    shutil.rmtree(wkdir)
