#!/usr/bin/python
#
# -- job name --
#
#$ -N /work/funashima/Si2_227
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

wkdir = '/work/funashima/Si2_227'
bindir = '/home/funashima/QE/espresso/bin'
inputfile = 'espresso_relax.in'
if os.path.isdir(wkdir):
    shutil.rmtree(wkdir)
os.mkdir(wkdir)

command = "{0}/pw.x < {1}".\
    format(bindir, inputfile)
proc = subprocess.call(command, shell=True)
if proc < 1:
    shutil.rmtree(wkdir)
