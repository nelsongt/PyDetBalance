from __future__ import division
import os
import math
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from scipy import integrate


## Some Constants ##

T_s = 6000 # K
T_c = 300 # K
k_B = 8.6173324E-05 # eV/K
planck_ev = 4.135667516E-15 # eV*s
subtend_angle_degrees = 0.26 # degrees

geometrical_factor = math.pi * math.sin(math.radians(subtend_angle_degrees)) ** 2

flux_coeff = 2 * geometrical_factor / (planck_ev**3 * sc.c**2)

## Integrate to find the flux ##

planck_rad_func = lambda E: E**3 / (math.exp(E/(k_B*T_s)) - 1)
flux, flux_err = integrate.quad(planck_rad_func,0,10) # Integrate the planck radiation function from 0 eV to 10 eV (~inf)


for x in xrange(1,20):
  y = (x/10)**2 / (math.exp((x/10)/(k_B*T_s)) - 1)
  print "Y: %g" %y


plt.plot(x, y, 'o', c = 'b')

print "Flux: %g" %flux