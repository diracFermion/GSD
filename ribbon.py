#!/usr/bin/python
import numpy as np
import hoomd
from hoomd import md, deprecated, data, group, init
import gsd.hoomd
import sys
import random

hoomd.context.initialize()
s = hoomd.init.read_gsd("../Sim_dump/init_strip.gsd")

for p in s.particles:
	if p.type == 'A':
		x, y, z = p.position
		z += random.uniform(-0.10,0.10)
		p.position = (x,y,z)

harmonic = md.bond.harmonic()
dih = md.dihedral.harmonic()

dih.dihedral_coeff.set('A', k=5.000, d=1, n=1)

harmonic.bond_coeff.set('A', k=800.000, r0=1.0)

hoomd.analyze.log(filename="../Sim_dump/observable.log", quantities=["temperature", "potential_energy","bond_harmonic_energy","kinetic_energy","dihedral_harmonic_energy"], period=5000, header_prefix="#", overwrite=True)

md.integrate.mode_standard(dt=0.0010)

group1 = hoomd.group.type(name='group1', type='A')

hoomd.dump.gsd(filename="../Sim_dump/trajectory.gsd", group=group.all(), period=5000, overwrite=True,static=[])

#print(s.particles[5])

md.integrate.nvt(group=group1,kT=1.0, tau=0.2)

#print(s.particles[5])

hoomd.run(10000000)

#print(s.particles[5])

#hoomd.run(10)
#print(s.particles[5])
#print(s.bonds[5])
