#!/usr/bin/python
#
# -- job name --
#
#$ -N Si2_221_1
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

wkdir = '/work/funashima/Si2_221_1'
bindir = '/home/funashima/QE/qe-6.4.1/bin'
inputfile = 'espresso_relax.in'
if os.path.isdir(wkdir):
    shutil.rmtree(wkdir)
os.mkdir(wkdir)

command = "{0}/pw.x < {1}".\
    format(bindir, inputfile)
proc = subprocess.call(command, shell=True)
if proc < 1:
    shutil.rmtree(wkdir)
