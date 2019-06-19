#!/usr/bin/env python3
#
#
#
from fractions import Fraction

class StrFraction(object):
    def __init__(self, n):
        self.n = n

    def add(self, b):
        self.n = str(Fraction(self.n) + Fraction(b))
        return self.n

    def multiply(self, b):
        self.n = str(Fraction(self.n) * Fraction(b))
        return self.n

    def show(self):
        return self.n


a = StrFraction('1/3')
print(a.show())
a.add('1/2')
print(a.show())
a.multiply('1/2')
print(a.show())
