#!/usr/bin/env python3
from Metis.SpaceGroup.ParseWyckoff import ParseWyckoff
from Metis.Structure.GenerateCrystal import GenerateAtomicPosition
from Metis.Structure.GetRandomLatticeConstant import GetRandomLatticeConstant
from Metis.SpaceGroup.GenerateSpaceGroup import GenerateSpaceGroup
import random

ispg = 'Fd-3m'
ichoice = 2
const_volume = 1000
max_coa_ratio = 4.0

random.seed(a=None, version=2)

space_group = GenerateSpaceGroup(ispg=ispg, ichoice=ichoice)
jspg = space_group.ispg # index of space group
crystal_system = space_group.crystal_system

lattice = GetRandomLatticeConstant(const_volume=const_volume,
                                   crystal_system=crystal_system)

a, b, c = lattice.get_lattice_length(max_coa_ratio=max_coa_ratio)
alpha, beta, gamma = lattice.get_lattice_angle()
print(a, b, c, alpha, beta, gamma)
lattice_params = {'a': a, 'b': b, 'c': c,
                  'alpha': alpha, 'beta': beta, 'gamma': gamma}
wyckoff =  ParseWyckoff()
wletters = ['d', 'h', 'i']
#for wlet in wletters:
wlet = 'h'
x = None
y = None
z = None
pos = wyckoff.get_atomic_position(ispg=jspg, ichoice=ichoice,
                                      wyckoff_letter = wlet)
print('** letter:{}'.format(wlet))
if any(['x' in x for x in pos]):
    x = random.random()
else:
    x = None
if any(['y' in x for x in pos]):
    y = random.random()
else:
    y = None
if any(['z' in x for x in pos]):
    z = random.random()
else:
    z = None
wyckoff_params = {'letter': wlet, 'variables':[x, y, z]}


obj = GenerateAtomicPosition(ispg=ispg, ichoice=ichoice,
                             wyckoff_params=wyckoff_params)
obj.show_info()

