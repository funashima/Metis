#!/usr/bin/env python3
#
# ref https://docs.python.org/ja/3/library/random.html
#
# generate lattice constant, a, b, c, alpha, beta and gamma automatically
# written by Hiroki Funashima in Kobe, 4 June 2019
#

import math
from Metis.Base.TspaceToolbox import TspaceToolbox
import random
import re


class GetRandomLatticeConstant(TspaceToolbox):
    def __init__(self, const_volume=None, crystal_system=None,
                 min_angle=60.0, max_angle=120.0):
        self.const_volume = const_volume
        self.crystal_system = crystal_system.lower()
        self.min_angle = min_angle
        self.max_angle = max_angle
        random.seed(a=None, version=2)

    def set_angle(self):
        self.alpha = 90.0
        self.beta = 90.0
        self.gamma = 90.0
        if re.search('^(cu|tet|ortho)', self.crystal_system):
            if re.search('^cu', self.crystal_system):
                self.freedom = 1
            elif re.search('^tet', self.crystal_system):
                self.freedom = 2
            else:
                self.freedom = 3
        elif re.search('^(rhom|hex|trig)', self.crystal_system):
            self.gamma = 120.0
            self.freedom = 2
        elif re.search('^mono', self.crystal_system):
            self.beta = self.random_deg_angle(min_angle=self.min_angle,
                                              max_angle=self.max_angle)
            self.freedom = 3
        elif re.search('^tric', self.crystal_system):
            self.alpha = self.random_deg_angle(min_angle=self.min_angle,
                                               max_angle=self.max_angle)
            self.beta = self.random_deg_angle(min_angle=self.min_angle,
                                              max_angle=self.max_angle)
            self.gamma = self.random_deg_angle(min_angle=self.min_angle,
                                               max_angle=self.max_angle)
            self.freedom = 3
        else:
            print('===== Error(set_angle) =====')
            print('unknown crystal system:{}'.format(self.crystal_system))
            exit()
        self.set_volume_coefficient()

    def set_volume_coefficient(self):
        """
          volume = a * b * c * self.vo
          this method, we get `self.vo'
          in detail, see notes
        """
        ca = self.deg2rad(self.gamma)
        cb = self.deg2rad(self.beta)
        cc = self.deg2rad(self.alpha)
        # ref http://gisaxs.com/index.php/Unit_cell
        try:
            self.vo = math.sqrt(1.0 + (2*ca*cb*cc) - (ca**2 + cb**2 + cc**2))
        except ValueError:  # in this case, self.vo is complex number.
            #
            # re-generate lattice angles alpha, beta and gamma
            #
            self.set_angle()

    def random_deg_angle(self, min_angle=60.0, max_angle=120.0):
        return random.uniform(min_angle, max_angle)

    def _get_lat_type2(self, mu, sigma, variate):
        a = -1.0
        while a < 0.0:
            if variate == 'normal':
                a = random.normalvariate(mu, sigma)
            elif variate == 'random':
                a = mu * random.random()
        b = self.const_volume / (a * a * self.vo)
        return [a, a, b]

    def _get_lat_type3(self, mu, sigma, variate):
        a = -1.0
        while a < 0.0:
            if variate == 'normal':
                a = random.normalvariate(mu, sigma)
            elif variate == 'random':
                a = mu * random.random()
        b = -1.0
        while b < 0.0:
            if variate == 'normal':
                b = random.normalvariate(mu, sigma)
            elif variate == 'random':
                b = mu * random.random()
        c = self.const_volume / (a * b * self.vo)
        return [a, b, c]

    def get_lattice_length(self,
                           max_coa_ratio=2.00,
                           variate='normal'):
        self.set_angle()
        coa_check = False
        while not coa_check:
            lat = self.random_length(variate)
            if max_coa_ratio > max(lat) / min(lat):
                coa_check = True
        return lat

    def random_length(self, variate):
        lat = [None] * 3
        r = self.const_volume ** (1/3)
        mu = r
        sigma = 1.0

        if self.freedom == 1:
            a = r
            lat = [a, a, a]
        elif self.freedom == 2:
            lat = self._get_lat_type2(mu, sigma, variate)
        elif self.freedom == 3:
            lat = self._get_lat_type3(mu, sigma, variate)
        else:
            print('===== Error(random_length) =====')
            print('freedom is incorrect. freedom = {}'.format(self.freedom))
            exit()
        return lat

    def get_lattice_angle(self):
        return [self.alpha, self.beta, self.gamma]

    def volume(self, a, b, c):
        return self.vo * a * b * c
