from __future__ import division
import os
import math
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from scipy import integrate


## Some Constants ##
T_s = 5760 # K
T_c = 300 # K
k_B = 8.6173324E-05 # eV/K
planck_ev = 4.135667516E-15 # eV*s
subtend_angle_degrees = 0.26 # degrees

sun_geo_factor = math.pi * math.sin(math.radians(subtend_angle_degrees)) ** 2
cell_geo_factor = math.pi

flux_const = 2 / (planck_ev**3 * sc.c**2)
#flux_const = 1 #cancels out for efficiency


## Required functions ##
def bose_einstein(E,T,qV):
  return 1 / (math.exp((E-qV)/(k_B*T)) - 1)

def planck_rad_func(E,T,qV):
  return E**3 * bose_einstein(E,T,qV)

def flux_func(E,T,qV):
  return E**2 * bose_einstein(E,T,qV)

def net_flux_func(E,qV):
  return sun_geo_factor*flux_func(E,T_s,0) - cell_geo_factor*flux_func(E,T_c,qV) + cell_geo_factor*flux_func(E,T_c,0)


#### MAIN ####

## Find the solar flux ##
y = np.zeros(50)
x = np.zeros(50)
for i in xrange(1,50):
  x[i] = i/10
  y[i] = flux_func(x[i],T_s,0) * flux_const * sun_geo_factor


plt.plot(x, y, 'o', c = 'b')
plt.show()


## Integrate to find the power from the sun ##
power_in, power_err = integrate.quad(planck_rad_func,0,10,args=(T_s,0)) # Integrate the planck radiation function from 0 eV to 10 eV (~inf)

power_in = power_in * flux_const * sun_geo_factor 

print "Power in: %g" %power_in


## Loop over voltages for Eg=1eV ##
y = np.zeros(50)
x = np.zeros(50)
z = np.zeros(50)
for i in xrange (1,50):
  x[i] = i/50
  ## Integrate to find current density for Eg=1ev ##
  J_0, J_err = integrate.quad(net_flux_func,1,10,args=(x[i]))
  y[i] = J_0 * flux_const
  z[i] = x[i]*y[i]
  
max = np.amax(z)/power_in
print "Max: %g" %max


## Loop over bandgaps ##
y = np.zeros(50)
x = np.zeros(50)
z = np.zeros(50)
E = np.zeros(40)
max = np.zeros(40)
for j in xrange (1,40):
  E[j] = j/10
  for i in xrange (1,50):
    x[i] = (i/50)*E[j]
    ## Integrate to find current density for Eg ##
    J_0, J_err = integrate.quad(net_flux_func,E[j],10,args=(x[i]))
    y[i] = J_0 * flux_const
    z[i] = x[i]*y[i]
  max[j] = np.amax(z)/power_in

plt.plot(E, max, 'o', c = 'b')
plt.show()