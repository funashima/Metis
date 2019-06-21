#!/usr/bin/env python3
import re

'''
   Parser Utility
     written by Hiroki Funashima in Kobe 2018-2019
'''


class ParseConfigfileUtil(object):
    def __init__(self):
        pass

    def _get_linebuf(self, line):
        linebuf = line.strip().replace(',', '')
        if '#' in linebuf:
            linebuf = linebuf.split('#')[0].strip()
        return linebuf

    def _parse_keyword_and_value(self, var_name, value, vtype=None):
        if vtype == 'int':
            exec('self.{0} = int({1})'.format(var_name, value))
        elif vtype == 'float':
            exec('self.{0} = float({1})'.
                 format(var_name,
                        value.lower().replace('d', 'e')))
        elif vtype == 'bool':
            if re.search('^(y|t)', value.lower()):
                exec('self.{} = True'.format(var_name, value))
            elif re.search('^(n|f)', value.lower()):
                exec('self.{} = False'.format(var_name, value))
        elif vtype == 'str':
            exec('self.{0} = "{1}"'.format(var_name, value))
        elif vtype == 'int_list':
            exec('self.{0} = list(map(lambda x:int(x), "{1}".split()))'.
                 format(var_name, value))
        elif vtype == 'float_list':
            exec('self.{0} = list(map(lambda x:float(x), "{1}".split()))'.
                 format(var_name, value))
        elif vtype == 'str_list':
            exec('self.{0} = list(map(lambda x:str(x), "{1}".split()))'.
                 format(var_name, value))
        else:
            exec('self.{0} = "{1}"'.format(var_name, value))

    def _get_key_and_value(self, linebuf):
        key, value = list(map(lambda x: x.strip(),
                              linebuf.split('=')[0:2]))
        key = key.lower()
        return [key, value]

    def _check_keyword_list(self, keyword_list, linebuf):
        key, value = self.get_key_and_value(linebuf)
        if value == '':
            return

        for (keyword, vtype) in keyword_list.items():
            if type(vtype) is list:
                if len(vtype) > 1:
                    val_name, val_type = vtype[:2]
                    if keyword == key:
                        self.parse_keyword_and_value(val_name, value, val_type)
                else:
                    if keyword == key:
                        self.parse_keyword_and_value(keyword, value, vtype[0])
            else:
                if keyword == key:
                    self.parse_keyword_and_value(keyword, value, vtype)

    def _check_data_region(self, region_name, linebuf):
        if '_' in linebuf:
            data = linebuf.split('_')
            key = data[0].lower().strip()
            value = data[1].lower().strip()
            if re.search('^begin', key) \
                    and re.search('^{}'.format(region_name), value):
                return True
            if re.search('^end', key) \
                    and re.search('^{}'.format(region_name), value):
                return False
            return None
        else:
            return None
