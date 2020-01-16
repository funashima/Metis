#!/usr/bin/python
#
# -- job name --
#
#$ -N Si2_ref
#
#$ -pe fillup 6
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

wkdir = '/work/funashima/ref'
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
