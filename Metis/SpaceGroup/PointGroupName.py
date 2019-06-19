#!/usr/bin/env python3
#
# class for PointGroup
#    written by HF in Kobe 28/07/2017
#    modified by HF in Kobe 19/06/2019
#
#    this program is based on Fortran version
#    written by HF 2016
#


class PointGroupName(object):
    def __init__(self, il=None, element_list=None, generator_list=None):
        self.element_list = element_list
        self.generator_list = generator_list
        self.crystal_system = None
        self.il = il

        if self.il > 0:
            self.__add_inv = 24
        else:
            self.__add_inv = 12

        c3axis, ic3axis = self.check_c3axes()
        if c3axis == 4:
            self.crystal_system = 1
            if self.check_cubic():
                self.name_cubic()
            else:
                self.catch_error(self.crystal_system)

        elif c3axis == 1:
            if any([self.check_c6(), self.check_ic6()]):
                self.crystal_system = 0
                if self.check_hexagonal():
                    self.name_hexagonal()
                else:
                    self.catch_error(self.crystal_system)
            else:
                self.crystal_system = -1
                if self.check_rhombohedral():
                    self.name_rhombohedral()
                else:
                    self.catch_error(self.crystal_system)

        elif c3axis == 0:
            c4axis, ic4axis = self.check_c4axes()
            if c4axis > 0 or ic4axis > 0:
                self.crystal_system = 2
                if self.check_tetragonal():
                    self.name_tetragonal()
                else:
                    self.catch_error(self.crystal_system)

            else:
                if self.check_orthorhombic():
                    self.crystal_system = 3
                    self.name_orthorhombic()
                else:
                    if self.check_monoclinic():
                        self.crystal_system = 4
                        self.name_monoclinic()
                    else:
                        if self.check_triclinic():
                            self.crystal_system = 5
                            self.name_triclinic()
                        else:
                            self.catch_error(9)
        else:
            self.catch_error(9)

    def catch_error(self, icode):
        print("========== Error(PointGroupName) ============")
        print(" undefined point group ")
        if icode == -1:
            print(' at rhombohedral')
        elif icode == 0:
            print(' at hexagonal')
        elif icode == 1:
            print(' at cubic')
        elif icode == 2:
            print(' at tetragonal')
        elif icode == 3:
            print(' at orthorhombic')
        elif icode == 4:
            print(' at monoclinic')
        elif icode == 5:
            print(' at trilinic')
        elif icode == 9:
            print(' at no crystal system')
        exit()

    def name_cubic(self):
        self.init_pg_name()
        if self.check_c4x():
            if self.check_inversion():
                self.pg_name = 'Oh'
            else:
                self.pg_name = 'O'
        else:
            if self.check_ic4x():
                self.pg_name = 'Td'
            else:
                if self.check_inversion():
                    self.pg_name = 'Th'
                else:
                    self.pg_name = 'T'

    def check_c4axes(self):
        c4axis = 0
        ic4axis = 0
        if self.check_c4x():
            c4axis += 1
        if self.check_c4y():
            c4axis += 1
        if self.check_c4z():
            c4axis += 1

        if self.check_ic4x():
            ic4axis += 1
        if self.check_ic4y():
            ic4axis += 1
        if self.check_ic4z():
            ic4axis += 1

        return [c4axis, ic4axis]

    def check_c2axes(self):
        if self.il > 0:
            c2axis, ic2axis = self.check_c2axes_Oh()
        else:
            c2axis, ic2axis = self.check_c2axes_D6h()
        return [c2axis, ic2axis]

    def check_c2axes_Oh(self):
        c2axis = self.check_c2_Oh()
        ic2axis = self.check_ic2_Oh()
        return [c2axis, ic2axis]

    def check_c2_Oh(self):
        c2axis = 0
        if self.check_c2x():
            c2axis += 1
        if self.check_c2y():
            c2axis += 1
        if self.check_c2z():
            c2axis += 1
        if self.check_c2a():
            c2axis += 1
        if self.check_c2b():
            c2axis += 1
        if self.check_c2c():
            c2axis += 1
        if self.check_c2d():
            c2axis += 1
        if self.check_c2e():
            c2axis += 1
        if self.check_c2f():
            c2axis += 1
        return c2axis

    def check_ic2_Oh(self):
        ic2axis = 0
        if self.check_ic2x():
            ic2axis += 1
        if self.check_ic2y():
            ic2axis += 1
        if self.check_ic2z():
            ic2axis += 1
        if self.check_ic2a():
            ic2axis += 1
        if self.check_ic2b():
            ic2axis += 1
        if self.check_ic2c():
            ic2axis += 1
        if self.check_ic2d():
            ic2axis += 1
        if self.check_ic2e():
            ic2axis += 1
        if self.check_ic2f():
            ic2axis += 1
        return ic2axis

    def check_c2axes_D6h(self):
        c2axis = self.check_c2_D6h()
        ic2axis = self.check_ic2_D6h()
        return [c2axis, ic2axis]

    def check_c2_D6h(self):
        c2axis = 0
        if self.check_c2():
            c2axis += 1
        if self.check_c211():
            c2axis += 1
        if self.check_c221():
            c2axis += 1
        if self.check_c231():
            c2axis += 1
        if self.check_c212():
            c2axis += 1
        if self.check_c222():
            c2axis += 1
        if self.check_c232():
            c2axis += 1
        return c2axis

    def check_ic2_D6h(self):
        ic2axis = 0
        if self.check_ic2():
            ic2axis += 1
        if self.check_ic211():
            ic2axis += 1
        if self.check_ic221():
            ic2axis += 1
        if self.check_ic231():
            ic2axis += 1
        if self.check_ic212():
            ic2axis += 1
        if self.check_ic222():
            ic2axis += 1
        if self.check_ic232():
            ic2axis += 1
        return ic2axis

    def check_c3axes(self):
        if self.il > 0:
            c3axis, ic3axis = self.check_c3_Oh()
        else:
            c3axis, ic3axis = self.check_c3_D6h()
        return [c3axis, ic3axis]

    def check_c3_Oh(self):
        c3axis = 0
        ic3axis = 0
        if self.check_c31():
            c3axis += 1
        if self.check_c32():
            c3axis += 1
        if self.check_c33():
            c3axis += 1
        if self.check_c34():
            c3axis += 1

        if self.check_ic31():
            ic3axis += 1
        if self.check_ic32():
            ic3axis += 1
        if self.check_ic33():
            ic3axis += 1
        if self.check_ic34():
            ic3axis += 1
        return [c3axis, ic3axis]

    def check_c3_D6h(self):
        c3axis = 0
        ic3axis = 0
        if self.check_c3():
            c3axis += 1
        if self.check_ic3():
            ic3axis += 1
        return [c3axis, ic3axis]

    def name_tetragonal(self):
        c2axis, ic2axis = self.check_c2axes()
        c4axis, ic4axis = self.check_c4axes()
        if self.check_inversion():
            if c2axis == 5:
                self.pg_name = 'D4h'
            elif c2axis == 1 and ic2axis == 1:
                self.pg_name = 'C4h'
        else:
            if c2axis == 5:
                self.pg_name = 'D4'
            elif c2axis == 1 and ic2axis == 4:
                self.pg_name = 'C4v'
            elif c2axis == 3 and ic2axis == 2 and ic4axis == 1:
                self.pg_name = 'D2d'
            elif c2axis == 1 and ic2axis != 4:
                if self.check_s4():
                    self.pg_name = 'S4'
                else:
                    self.pg_name = 'C4'

    def name_orthorhombic(self):
        self.init_pg_name()
        if self.check_inversion():
            self.pg_name = 'D2h'
        else:
            c2axis, ic2axis = self.check_c2axes()
            if c2axis > 2 and ic2axis == 0:
                self.pg_name = 'D2'
            elif c2axis > 0 and ic2axis > 1:
                self.pg_name = 'C2v'
            else:
                print("==== Error(PointGroupName) ====")
                print(" unknown point group(orthorhobmic)")
                print(" logic error")
                exit()

    def name_monoclinic(self):
        self.init_pg_name()
        if self.check_inversion():
            self.pg_name = 'C2h'
        else:
            c2axis, ic2axis = self.check_c2axes()
            if c2axis > 0 and ic2axis == 0:
                self.pg_name = 'C2'
            elif c2axis == 0 and ic2axis > 0:
                self.pg_name = 'Cs'
            else:
                print("==== Error(PointGroupName) ====")
                print(" unknown point group(monoclinic)")
                print(" logic error")
                print(" IL = {0}".format(self.il))
                print(" generator:", end='')
                print(self.generator_list)
                exit()

    def name_triclinic(self):
        self.init_pg_name()
        if self.check_inversion():
            self.pg_name = 'Ci'
        else:
            self.pg_name = 'C1'

    def name_hexagonal(self):
        self.init_pg_name()
        if self.check_inversion():
            if self.check_c211():
                self.pg_name = 'D6h'
            else:
                self.pg_name = 'C6h'
        else:
            if self.check_ic6():
                if self.check_ic2():
                    if any([self.check_c211(), self.check_c212()]):
                        self.pg_name = 'D3h'
                    else:
                        self.pg_name = 'C3h'
                else:
                    if self.check_ic211():
                        self.pg_name = 'C6v'
                    elif self.check_ic2():
                        self.pg_name = 'C6h'
                    else:
                        print("==== Error(PointGroupName) ====")
                        print(" unknown point group(hexagonal)")
                        print(" logic error")
                        exit()
            else:
                if self.check_c211():
                    self.pg_name = 'D6'
                else:
                    if any([self.check_ic211(), self.check_ic221()]):
                        self.pg_name = 'C6v'
                    else:
                        self.pg_name = 'C6'

    def name_rhombohedral(self):
        self.init_pg_name()
        c2axis, ic2axis = self.check_c2axes()
        c3axis, ic3axis = self.check_c3axes()
        if c3axis == 0:
            print("==== Error(PointGroupName) ====")
            print(" unknown point group(rhombohedral)")
            print(" logic error")
            exit()
        if self.check_inversion():
            if c2axis == 3:
                self.pg_name = 'D3d'
            elif c2axis == 0:
                self.pg_name = 'C3i'
            else:
                print("==== Error(PointGroupName) ====")
                print(" unknown point group(rhombohedral)")
                print(" logic error")
                exit()

        else:
            if ic2axis == 3:
                self.pg_name = 'C3v'
            elif c2axis == 3:
                self.pg_name = 'D3'
            elif c2axis == 0 and ic2axis == 0:
                self.pg_name = 'C3'
            else:
                print("==== Error(PointGroupName) ====")
                print(" unknown point group(rhombohedral)")
                print(" logic error")
                exit()

    def init_pg_name(self):
        self.pg_name = ''
        return self

    def check_ig(self, icode):
        if icode in self.element_list:
            return True
        else:
            return False

    def check_ng(self):
        return len(self.element_list)

    def check_inversion(self):
        if self.il > 0:
            icode = 25
        else:
            icode = 13
        return self.check_ig(icode)

    #
    # for cubic system
    #
    def check_cubic(self):
        return self.check_c31()

    #
    # for tetragonal system
    #

    def check_tetragonal(self):
        if self.check_c31():
            return False
        else:
            if any([self.check_c4z(), self.check_ic4z()]):
                return True
            else:
                return False

    #
    # for Orthorhobmic System
    #

    def check_orthorhombic(self):
        if not any([self.check_cubic(), self.check_tetragonal()]):
            if self.check_c2z():
                if any([self.check_c2x(),
                        self.check_ic2x(),
                        self.check_ic2a()]):
                    return True
                else:
                    return False
            else:
                if all([self.check_c2y(), self.check_ic2x()]):
                    return True
                else:
                    return False
        else:
            return False

    def check_monoclinic(self):
        if any([self.check_cubic(),
                self.check_tetragonal(),
                self.check_orthorhombic()]):
            return False
        else:
            if self.check_c2y() is not self.check_c2z():
                return True
            elif self.check_ic2y() is not self.check_ic2z():
                return True
            else:
                return False

    def check_triclinic(self):
        if any([self.check_cubic(),
                self.check_tetragonal(),
                self.check_orthorhombic(),
                self.check_monoclinic()]):
            return False
        elif any([self.check_c2y(),
                  self.check_c2z(),
                  self.check_ic2y(),
                  self.check_ic2z()]):
            return False
        else:
            return True

    def check_hexagonal(self):
        #
        # this hexagonal is not included D3d
        #
        if any([self.check_c6(), self.check_ic6()]):
            return True
        else:
            return False

    def check_rhombohedral(self):
        if any([self.check_c6(), self.check_ic6()]):
            return False
        else:
            if self.check_c3():
                return True
            else:
                return False

    def check_c2x(self):
        icode = 2
        return self.check_Oh_subgp(icode)

    def check_c2y(self):
        icode = 3
        return self.check_Oh_subgp(icode)

    def check_c2z(self):
        icode = 4
        return self.check_Oh_subgp(icode)

    def check_c31(self):
        icode = 5
        return self.check_Oh_subgp(icode)

    def check_c32(self):
        icode = 6
        return self.check_Oh_subgp(icode)

    def check_c33(self):
        icode = 7
        return self.check_Oh_subgp(icode)

    def check_c34(self):
        icode = 8
        return self.check_Oh_subgp(icode)

    # (9-12) :    c31- -> c34-

    def check_c2a(self):
        icode = 13
        return self.check_Oh_subgp(icode)

    def check_c2b(self):
        icode = 14
        return self.check_Oh_subgp(icode)

    def check_c2c(self):
        icode = 15
        return self.check_Oh_subgp(icode)

    def check_c2d(self):
        icode = 16
        return self.check_Oh_subgp(icode)

    def check_c2e(self):
        icode = 17
        return self.check_Oh_subgp(icode)

    def check_c2f(self):
        icode = 18
        return self.check_Oh_subgp(icode)

    def check_c4x(self):
        icode = 19
        return self.check_Oh_subgp(icode)

    def check_c4y(self):
        icode = 20
        return self.check_Oh_subgp(icode)

    def check_c4z(self):
        icode = 21
        return self.check_Oh_subgp(icode)

    # (22-24): c4x- -> c4z-

    def check_ic2x(self):
        icode = 2 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2y(self):
        icode = 3 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2z(self):
        icode = 4 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic31(self):
        icode = 5 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic32(self):
        icode = 6 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic33(self):
        icode = 7 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic34(self):
        icode = 8 + self.__add_inv
        return self.check_Oh_subgp(icode)

    # (9-12) + 24:  ic31- -> ic34-

    def check_ic2a(self):
        icode = 13 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2b(self):
        icode = 14 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2c(self):
        icode = 15 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2d(self):
        icode = 16 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2e(self):
        icode = 17 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic2f(self):
        icode = 18 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic4x(self):
        icode = 19 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic4y(self):
        icode = 20 + self.__add_inv
        return self.check_Oh_subgp(icode)

    def check_ic4z(self):
        icode = 21 + self.__add_inv
        return self.check_Oh_subgp(icode)

    # (22 - 24) + 24 :  ic4x- -> ic4z-

    #
    # function for D6h subgroup
    #

    def check_c6(self):
        icode = 2
        return self.check_D6h_subgp(icode)

    def check_c3(self):
        icode = 3
        return self.check_D6h_subgp(icode)

    def check_c2(self):
        icode = 4
        return self.check_D6h_subgp(icode)

    def check_c211(self):
        icode = 7
        return self.check_D6h_subgp(icode)

    def check_c221(self):
        icode = 8
        return self.check_D6h_subgp(icode)

    def check_c231(self):
        icode = 9
        return self.check_D6h_subgp(icode)

    def check_c212(self):
        icode = 10
        return self.check_D6h_subgp(icode)

    def check_c222(self):
        icode = 11
        return self.check_D6h_subgp(icode)

    def check_c232(self):
        icode = 12
        return self.check_D6h_subgp(icode)

    def check_ic6(self):
        icode = 2 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic3(self):
        icode = 3 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic2(self):
        icode = 4 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic211(self):
        icode = 7 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic221(self):
        icode = 8 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic231(self):
        icode = 9 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic212(self):
        icode = 10 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic222(self):
        icode = 11 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_ic232(self):
        icode = 12 + self.__add_inv
        return self.check_D6h_subgp(icode)

    def check_s4(self):
        check = False
        if self.check_c2x() and self.check_ic4x():
            check = True
        if self.check_c2y() and self.check_ic4y():
            check = True
        if self.check_c2z() and self.check_ic4z():
            check = True
        return check

    def check_Oh_subgp(self, icode):
        if self.il < 1:
            return False
        else:
            return self.check_ig(icode)

    def check_D6h_subgp(self, icode):
        if self.il > 0:
            return False
        else:
            return self.check_ig(icode)
