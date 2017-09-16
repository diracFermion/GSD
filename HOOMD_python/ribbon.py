#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
import sys
import random

hoomd.context.initialize()
s = hoomd.init.read_gsd("../../Sim_dump_ribbon/init_strip.gsd")

for p in s.particles:
	if p.type == 'A':
		x, y, z = p.position
		z += random.uniform(-0.10,0.10)
		p.position = (x,y,z)

harmonic = md.bond.harmonic()
dih = md.dihedral.harmonic()

dih.dihedral_coeff.set('A', k=4.000, d=1, n=1)
dih.dihedral_coeff.set('D', k=4.000, d=1, n=1)
dih.dihedral_coeff.set('E', k=4.000, d=1, n=1)

harmonic.bond_coeff.set('A', k=800.000, r0=1.0)
harmonic.bond_coeff.set('D', k=800.000, r0=1.0)
harmonic.bond_coeff.set('E', k=800.000, r0=1.0)

hoomd.analyze.log(filename="../../Sim_dump_ribbon/Observable.log", quantities=["temperature", "potential_energy","bond_harmonic_energy","kinetic_energy","dihedral_harmonic_energy"], period=5000, header_prefix="#", overwrite=True)

md.integrate.mode_standard(dt=0.0010)

group1 = hoomd.group.type(name='group1', type='A')
group2 = hoomd.group.type(name='group2', type='D')
group3 = hoomd.group.type(name='group3', type='E')
group12 = hoomd.group.union(name='group12',a=group1,b=group2)
group123 = hoomd.group.union(name='group123',a=group12,b=group3)

md.constrain.oneD(group=group2, constraint_vector=[1,0,0])

hoomd.dump.gsd(filename="../../Sim_dump_ribbon/trajectory1.gsd", group=group.all(), period=5000, overwrite=True)


md.integrate.nvt(group=group123,kT=1.0, tau=0.2)

hoomd.run(1e6)

