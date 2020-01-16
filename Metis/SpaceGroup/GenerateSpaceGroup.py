#!/usr/bin/env python3
from Metis.Base.TspaceToolbox import TspaceToolbox
from Metis.SpaceGroup.ParseGenerator import ParseGenerator
from Metis.SpaceGroup.PointGroupName import PointGroupName
from fractions import Fraction
import sys


class GenerateSpaceGroup(TspaceToolbox):
    def __init__(self, ispg, ichoice=1):
        Gen_obj = ParseGenerator()
        self.generator = Gen_obj.get_generator(ispg, ichoice)
        self.il = Gen_obj.get_il(ispg, ichoice)
        self.schname = Gen_obj.get_schname(ispg, ichoice)
        self.hmname = Gen_obj.get_hmname(ispg, ichoice)
        self.ispg = Gen_obj.get_ispg(ispg, ichoice)
        self.ichoice = ichoice
        self.define_rotation_operator()
        self.initialize_operator()
        self.generate_group_elements()
        self.set_point_group_info()

    def define_rotation_operator(self):
        if self.il in [0, -1]:
            self.define_rotation_operator_D6h()
            self.define_operator_name_D6h()
            return self
        if self.il in [1, 2, 3, 4]:
            self.define_rotation_operator_Oh()
            self.define_operator_name_Oh()
            return self
        else:
            print('==== error in define_rotation_operator ==== ')
            print('il:{} is incorrect'.format(self.il))
            exit()

    def get_crystal_system(self, crystal_system_code):
        if crystal_system_code == -1:
            return 'rhombohedral'
        elif crystal_system_code == 0:
            return 'hexagonal'
        elif crystal_system_code == 1:
            return 'cubic'
        elif crystal_system_code == 2:
            return 'tetragonal'
        elif crystal_system_code == 3:
            return 'orthorhombic'
        elif crystal_system_code == 4:
            return 'monoclinic'
        elif crystal_system_code == 5:
            return 'triclinic'
        else:
            print("====== Error(PointGroupName) ======")
            print(" unknown crystal system")
            exit()

    def get_bravais_lattice_name(self):
        if self.il <= 0:
            prefix = ''
        elif self.il == 1:
            prefix = 'simple'
        elif self.il == 2:
            prefix = 'face centered'
        elif self.il == 3:
            prefix = 'body centered'
        elif self.il == 4:
            prefix = 'bace centered'
        else:
            print('unknown lattice type')
            print('il = {0}'.format(self.il))
            exit()
        return prefix + ' ' + self.crystal_system + ' lattice'

    def set_point_group_info(self):
        generator_list = []
        for generator in self.generator:
            generator_list.append(generator['rotation'])
        pg_obj = PointGroupName(il=self.il,
                                element_list=self.point_group,
                                generator_list=generator_list)
        crystal_system_code = pg_obj.crystal_system
        self.crystal_system = self.get_crystal_system(crystal_system_code)
        self.bravais_lattice_name = self.get_bravais_lattice_name()
        self.point_group_name = pg_obj.pg_name

    def initialize_operator(self):
        if self.schname != 'C11':
            identical_op = {'rotation': 1, 'translation': [0, 1, 0, 1, 0, 1]}
            self._wk_group_elements = [identical_op]
            self.point_group = [1]
        else:
            self._wk_group_elements = []
            self.point_group = []
        for generator in self.generator:
            self._wk_group_elements.append(generator)
            self.point_group.append(generator['rotation'])

    def sort_wk_group_elements(self):
        self.group_elements = []
        if self.il in [0, -1]:
            max_gen = 24
        else:
            max_gen = 48
        for i in range(1, max_gen + 1):
            for element in self._wk_group_elements:
                op_name = self.operator_name_list[i-1]
                if element['rotation'] == i:
                    char_rep = self.tspcode2char(element['rotation'])
                    mtx = self.convert_rotmatrix(char_rep)
                    new_element = {'rotation': element['rotation'],
                                   'translation': element['translation'],
                                   'rot_matrix': mtx,
                                   'operator_name': op_name}
                    self.group_elements.append(new_element)

    def generate_group_elements(self):
        check_full = False
        while not check_full:
            ngen = len(self._wk_group_elements)
            get_new_element = False
            for i in range(ngen):
                for j in range(ngen):
                    gen1 = self._wk_group_elements[i]
                    gen2 = self._wk_group_elements[j]
                    new_gen = self.multiply_generate_operators(gen1, gen2)
                    new_pg = new_gen['rotation']
                    if new_pg not in self.point_group:
                        self._wk_group_elements.append(new_gen)
                        self.point_group.append(new_pg)
                        get_new_element = True
                        break
            if not get_new_element:
                check_full = True
        self.sort_wk_group_elements()
        self.modified_tvec()

    def modified_tvec(self):
        for (i, element) in enumerate(self.group_elements):
            tvec = element['translation']
            m_tvec = self.check_tvec(tvec)
            self.group_elements[i]['translation'] = m_tvec
        return self

    def check_tvec(self, tvec):
        if self.il in [0, 1]:
            return tvec
        u = self.tsp2frac(tvec)
        if self.il == 2:
            u = self.check_tvec_fc(u)
        elif self.il == 3:
            u = self.check_tvec_bc(u)
        elif self.il == 4:
            u = self.check_tvec_cc(u)
        elif self.il == -1:
            u = self.check_tvec_rh(u)
        return self.frac2tsp(u)

    def check_tvec_fc(self, u):
        half = Fraction('1/2')
        if all([x >= half for x in u[:2]]):
            for i in range(2):
                u[i] -= half
            return u
        if all([x >= half for x in u[1:3]]):
            for i in range(1, 3):
                u[i] -= half
            return u
        if all([x >= half for x in [u[0], u[2]]]):
            for i in [0, 2]:
                u[i] -= half
            return u
        return u

    def check_tvec_bc(self, u):
        half = Fraction('1/2')
        if all([x >= half for x in u]):
            for i in range(3):
                u[i] -= half
        return u

    def check_tvec_cc(self, u):
        half = Fraction('1/2')
        if all([x >= half for x in u[:2]]):
            for i in range(2):
                u[i] -= half
        return u

    def check_tvec_rh(self, u):
        tr1 = Fraction('1/3')
        tr2 = 2 * tr1
        if u[0] >= tr2 and all([x >= tr1 for x in u[1:3]]):
            u[0] -= tr2
            for i in [1, 2]:
                u[i] -= tr1
            return u
        if u[0] >= tr1 and all([x >= tr2 for x in u[1:3]]):
            u[0] -= tr1
            for i in [1, 2]:
                u[i] -= tr2
            return u
        return u

    def multiply_generate_operators(self, gen1, gen2):
        #
        # In Seitz natation
        #  (r2|t2) (r1|t1) = (r2r1| r2t1 + t2)
        #
        #
        rot1 = self.convert_rotmatrix(self.tspcode2char(gen1['rotation']))
        rot2 = self.convert_rotmatrix(self.tspcode2char(gen2['rotation']))
        tvec1 = gen1['translation']
        tvec2 = gen2['translation']
        #
        # get r2 * r1
        #
        new_rotation = self.multiply_matrices(rot2, rot1)
        new_rotation = self.char2tspcode(new_rotation)

        #
        # get `r2t1 + t2'
        #
        new_tvec = self.rotate_tvec(tvec1, tvec2, rot2)
        return {'rotation': new_rotation, 'translation': new_tvec}

    def rotate_tvec(self, tvec1, tvec2, rot):
        #
        # In Seitz natation
        #  (r2|t2) (r1|t1) = (r2r1| r2t1 + t2)
        #
        # in this method, you can get `r2t1 + t2'
        #
        #
        tvec = [0 for x in range(3)]
        tt1 = self.tsp2frac(tvec1)  # convert fractional notation
        tt2 = self.tsp2frac(tvec2)  # convert fractional notation
        for i in range(3):
            for j in range(3):
                tvec[i] += rot[i][j] * tt1[j]
            tvec[i] += tt2[i]
        for i in range(3):
            if tvec[i] < 0:
                tvec[i] = 1 + tvec[i]
            if tvec[i] > 1:
                tvec[i] = -1 + tvec[i]
        return self.frac2tsp(tvec)

    def tspcode2char(self, tspcode):
        return self.operator_list[tspcode - 1]

    def char2tspcode(self, matrix):
        operator = [None] * 3
        for i in range(3):
            operator[i] = self.to_charcode(matrix[i])
        tspcode = self.operator_list.index(operator)
        if tspcode is not None:
            return tspcode + 1
        else:
            print('===== Error =====')
            print('unknown operator in matrix2tspcode')
            exit()

    def define_rotation_operator_Oh(self):
        self.operator_list = [
            ['x', 'y', 'z'],
            ['x', '-y', '-z'],
            ['-x', 'y', '-z'],
            ['-x', '-y', 'z'],
            ['z', 'x', 'y'],
            ['-z', 'x', '-y'],
            ['-z', '-x', 'y'],
            ['z', '-x', '-y'],
            ['y', 'z', 'x'],
            ['y', '-z', '-x'],
            ['-y', 'z', '-x'],
            ['-y', '-z', 'x'],
            ['y', 'x', '-z'],
            ['-y', '-x', '-z'],
            ['z', '-y', 'x'],
            ['-x', 'z', 'y'],
            ['-z', '-y', '-x'],
            ['-x', '-z', '-y'],
            ['x', '-z', 'y'],
            ['z', 'y', '-x'],
            ['-y', 'x', 'z'],
            ['x', 'z', '-y'],
            ['-z', 'y', 'x'],
            ['y', '-x', 'z'],
            ['-x', '-y', '-z'],
            ['-x', 'y', 'z'],
            ['x', '-y', 'z'],
            ['x', 'y', '-z'],
            ['-z', '-x', '-y'],
            ['z', '-x', 'y'],
            ['z', 'x', '-y'],
            ['-z', 'x', 'y'],
            ['-y', '-z', '-x'],
            ['-y', 'z', 'x'],
            ['y', '-z', 'x'],
            ['y', 'z', '-x'],
            ['-y', '-x', 'z'],
            ['y', 'x', 'z'],
            ['-z', 'y', '-x'],
            ['x', '-z', '-y'],
            ['z', 'y', 'x'],
            ['x', 'z', 'y'],
            ['-x', 'z', '-y'],
            ['-z', '-y', 'x'],
            ['y', '-x', '-z'],
            ['-x', '-z', 'y'],
            ['z', '-y', '-x'],
            ['-y', 'x', '-z']
            ]

    def define_operator_name_Oh(self):
        self.operator_name_list = [
            'E',
            'C2X',  'C2Y',  'C2Z',
            'C31+', 'C32+', 'C33+', 'C34+',
            'C31-', 'C32-', 'C33-', 'C34-',
            'C2A',  'C2B',  'C2C',  'C2D',  'C2E',  'C2F',
            'C4X+', 'C4Y+', 'C4Z+',
            'C4X-', 'C4Y-', 'C4Z-',
            'IE',
            'IC2X',  'IC2Y',  'IC2Z',
            'IC31+', 'IC32+', 'IC33+', 'IC34+',
            'IC31-', 'IC32-', 'IC33-', 'IC34-',
            'IC2A',  'IC2B',  'IC2C',  'IC2D',  'IC2E',  'IC2F',
            'IC4X+', 'IC4Y+', 'IC4Z+',
            'IC4X-', 'IC4Y-', 'IC4Z-'
            ]

    def define_rotation_operator_D6h(self):
        self.operator_list = [
            ['x', 'y', 'z'],
            ['w', 'x', 'z'],
            ['-y', 'w', 'z'],
            ['-x', '-y', 'z'],
            ['-w', '-x', 'z'],
            ['y', '-w', 'z'],
            ['-w', 'y', '-z'],
            ['x', 'w', '-z'],
            ['-y', '-x', '-z'],
            ['w', '-y', '-z'],
            ['-x', '-w', '-z'],
            ['y', 'x', '-z'],
            ['-x', '-y', '-z'],
            ['-w', '-x', '-z'],
            ['y', '-w', '-z'],
            ['x', 'y', '-z'],
            ['w', 'x', '-z'],
            ['-y', 'w', '-z'],
            ['w', '-y', 'z'],
            ['-x', '-w', 'z'],
            ['y', 'x', 'z'],
            ['-w', 'y', 'z'],
            ['x', 'w', 'z'],
            ['-y', '-x', 'z']
            ]

    def define_operator_name_D6h(self):
        self.operator_name_list = [
            'E',
            'C6+', 'C3+',
            'C2',
            'C6-', 'C3-',
            'C211', 'C221', 'C231',
            'C212', 'C222', 'C232',
            'IE',
            'IC6+', 'IC3+',
            'IC2',
            'IC6-', 'IC3-',
            'IC211', 'IC221', 'IC231',
            'IC212', 'IC222', 'IC232',
            ]

    def convert_rotmatrix(self, code):
        matrix = []
        for x in code:
            matrix.append(self.to_matrix_representation(x))
        return matrix

    def display_group_table(self, filename=None):
        def sub_display_table(fout):
            fout.write('Group Table\n')
            ngen = len(self.group_elements)
            for i in range(ngen):
                for j in range(ngen):
                    gen1 = self.group_elements[i]
                    gen2 = self.group_elements[j]
                    new_gen = self.multiply_generate_operators(gen2, gen1)
                    fout.write(' {:>2d}'.format(new_gen['rotation']))
                fout.write('\n')
            fout.write('\n')

        if filename is None:
            sub_display_table(sys.stdout)
        else:
            with open(filename, 'a') as fout:
                sub_display_table(fout)
        return self

    def display_group_elements(self, filename=None):
        def sub_display_group_elements(fout):
            fout.write('Group Elements\n')
            for (i, element) in enumerate(self.group_elements):
                fout.write(' {:>2d}'.format(i+1))
                tspcode = element['rotation']
                fout.write('  {:>2d}'.format(tspcode))
                fout.write('  {0:6s}'.format(element['operator_name']))
                for x in self.tspcode2char(tspcode):
                    fout.write(' {:>2s}'.format(x.upper()))
                fout.write('  ')
                for i in range(3):
                    fout.write(' {0}/{1}'.
                               format(element['translation'][2*i],
                                      element['translation'][2*i+1]))
                    if i == 2:
                        fout.write('\n')
                    else:
                        fout.write(' ')
            fout.write('\n')
        if filename is None:
            sub_display_group_elements(sys.stdout)
        else:
            with open(filename, 'a') as fout:
                sub_display_group_elements(fout)
        return self
