#!/usr/bin/env python3
#
# convert pwscf inputfile to space group
#   written by Hiroki Funashima in Kobe, 30 May 2019
#   modified by Hiroki Funashima in Kobe, 3 June 2019
#

import os


class QECheckFileType(object):
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename
        else:
            print('===== Error(QECheckFileType) =====')
            print('file:{} is not found.'.format(filename))
            exit()
        self.main()

    def _get_linebuf(self, line):
        linebuf = line.strip().lower()
        for comment in ['#', '!']:
            if comment in linebuf:
                linebuf = linebuf.split(comment)[0].strip()
        return linebuf

    def main(self):
        self.filetype = 'optimize'
        for line in open(self.filename, 'r'):
            linebuf = self._get_linebuf(line)
            if 'begin' in linebuf:
                if 'final' in linebuf:
                    if 'coordinates' in linebuf:
                        self.filetype = 'optimize'
                        return
            if '&system' in linebuf:
                self.filetype = 'inputfile'
                return
        return 
