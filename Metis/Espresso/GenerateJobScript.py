#!/usr/bin/env python3
#
# generate job control script
#    written by Hiroki Funashima, 25 June, 2019
#
import os


class ParseConfig(object):
    def __init__(self, configfile):
        if os.path.isfile(configfile):
            self.configfile = configfile
        else:
            print('file:{} is not found.'.format(configfile))
            exit()
        self.main()

    def get_linebuf(self, line):
        linebuf = line.strip()
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def get_key_and_value(self, linebuf):
        try:
            key, value = [x.strip() for x in linebuf.split('=')]
        except:
            print(linebuf)
        key = key.lower()
        return [key, value]

    def main(self):
        self.nnode = 1
        for line in open(self.configfile, 'r'):
            linebuf = self.get_linebuf(line)
            if linebuf == '':
                continue
            if '=' not in linebuf:
                continue
            if ';' in linebuf:
                continue
            key, value = self.get_key_and_value(linebuf)
            if value == '':
                continue
            if key == 'bindir':
                self.bindir = value
            elif key == 'nnode':
                self.nnode = int(value)


class GenerateJobScript(object):
    def __init__(self, configfile=None,
                 wkdir=None,
                 calc_dir=None,
                 inputfile=None):
        self.configure = ParseConfig(configfile)
        self.wkdir = wkdir
        self.calc_dir = calc_dir
        self.inputfile = inputfile
        self.scriptfile = os.path.join(calc_dir, 'job.py')
        self.main()

    def main(self):
        self.write_header()
        self.write_imports()
        self.write_main()

    def write_header(self):
        nnode = self.configure.nnode
        with open(self.scriptfile, 'w') as fout:
            fout.write('#!/usr/bin/python\n')
            fout.write('#\n')
            fout.write('# -- job name --\n')
            fout.write('#\n')
            fout.write('#$ -N {}\n'.format(self.wkdir))
            fout.write('#\n')
            fout.write('#$ -pe fillup {}\n'.format(nnode))
            fout.write('#\n')
            fout.write('#$ -V\n')
            fout.write('#$ -S /usr/bin/python\n')
            fout.write('#$ -cwd\n')
            fout.write('#$ -e espresso_relax.err\n')
            fout.write('#$ -o espresso_relax.out\n')
            fout.write('#\n')
            fout.write('\n')
        return self

    def write_imports(self):
        with open(self.scriptfile, 'a') as fout:
            fout.write('import os\n')
            fout.write('import shutil\n')
            fout.write('import subprocess\n')
            fout.write('\n')
        return self

    def write_main(self):
        bindir = self.configure.bindir
        indent = '    '
        with open(self.scriptfile, 'a') as fout:
            fout.write("wkdir = '{}'\n".format(self.wkdir))
            fout.write("bindir = '{}'\n".format(bindir))
            fout.write("inputfile = '{}'\n".format(self.inputfile))
            fout.write('if os.path.isdir(wkdir):\n')
            fout.write('{}shutil.rmtree(wkdir)\n'.format(indent))
            fout.write('os.mkdir(wkdir)\n'.format(indent))
            fout.write('\n')
            fout.write('command = "{0}/pw.x < {1}".\\\n')
            fout.write('    format(bindir, inputfile)\n')
            fout.write('proc = subprocess.call(command, shell=True)\n')
            fout.write('if proc < 1:\n')
            fout.write('{}shutil.rmtree(wkdir)\n'.format(indent))
        return self
