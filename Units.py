#!/usr/bin/env python
#
# Useful physical units and constants
#
# Created by: dadwyer@lbl.gov 2013/03/06
# Last Modified: dadwyer@lbl.gov 2013/03/15
# 
#

# Units
# energy
eV = 1.0
MeV = 1.0e6 * eV
GeV = 1.0e9 * eV
# length
m = 1.0
fm = 1.0e-15 * m
cm = 1.0e-2 * m
km = 1.0e3 * m
# time
s = 1.0
year = 365*24*60*60 * s
# mass
g = 1.0
kg = 1e3 * g
ton = 1e6 * g
Mton = 1e6 * ton
# solid angle
steradian = 1

# Useful constants
hbarc = 197.326968 * MeV * fm
amu = 1.66053886e-27 * kg
earthRadius = 6378.140 * km
Gf = 1.16637e-5 * GeV**-2 * hbarc**3
rho_e = (4.5 * g / cm**3) * (0.5 / amu)
