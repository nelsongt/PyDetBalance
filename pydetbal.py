from __future__ import division
import os
import math
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize


## Some Constants ##
T_s = 5760 # K
T_c = 300 # K
k_B = 8.6173324E-05 # eV/K
planck_ev = 4.135667516E-15 # eV*s
subtend_angle_degrees = 0.26 # degrees
conc_factor = 1 # suns

sun_geo_factor = math.pi * math.sin(math.radians(subtend_angle_degrees)) ** 2
cell_geo_factor = math.pi

#flux_const = 2 / (planck_ev**3 * sc.c**2)
flux_const = 1 #cancels out for efficiency


## Required functions ##
def bose_einstein(E,T,qV):
  return 1 / (math.exp((E-qV)/(k_B*T)) - 1)

def planck_rad_func(E,T,qV):
  return E**3 * bose_einstein(E,T,qV)

def flux_func(E,T,qV):
  return E**2 * bose_einstein(E,T,qV)

def net_flux_func(E,qV):
  return conc_factor*sun_geo_factor*flux_func(E,T_s,0) - cell_geo_factor*flux_func(E,T_c,qV) + (cell_geo_factor - conc_factor*sun_geo_factor)*flux_func(E,T_c,0)

def current_from_flux(qV,E1,E2):
  J_0, J_err = integrate.quad(net_flux_func,E1,E2,args=(qV))
  return J_0

def net_current_IBSC(inter_fermi,qV,Eiv,Eic,Eg):
  return abs(current_from_flux(inter_fermi,Eiv,Eic) - current_from_flux(qV-inter_fermi,Eic,Eg))


#### MAIN ####

## Find the solar flux ##
#y = np.zeros(50)
#x = np.zeros(50)
#for i in xrange(1,50):
#  x[i] = i/10
#  y[i] = flux_func(x[i],T_s,0) * flux_const * sun_geo_factor

#plt.plot(x, y, 'o', c = 'b')
#plt.show()


## Integrate to find the power from the sun ##
power_in, power_err = integrate.quad(planck_rad_func,0,10,args=(T_s,0)) # Integrate the planck radiation function from 0 eV to 10 eV (~inf)

power_in = power_in * flux_const * sun_geo_factor * conc_factor

#print "Power in: %g" %power_in



## Loop over bandgaps ##
## This got complicated quickly and I should rewrite it to be more readable ##
V = np.zeros(50)
J = np.zeros(50)
P = np.zeros(50)
Vmax = np.zeros(40)
Eiv = np.zeros(40)
Eic = np.zeros(40)
Eg = np.zeros(40)
max_eff_1 = np.zeros(40)
max_Eic = np.zeros(40)
max_eff_2 = np.zeros(40)
for j in xrange (10,20):
  Eiv[j] = j/15 # Choose intermediate band energy first
  for i in xrange(5,25):
    Eg[i] = 2*Eiv[j] + i/20 # Choose bulk bandgap next, should be over twice the Ei value
    Eic[i] = Eg[i] - Eiv[j]
    for k in xrange(1,50):
      V[k] = Eiv[j] + k/50 * Eic[i] # Constrain voltages
      inter_fermi = optimize.fminbound(net_current_IBSC,V[k]-Eic[i],Eiv[j],args=(V[k],Eiv[j],Eic[i],Eg[i]))
      J[k] = flux_const * (current_from_flux(V[k],Eg[i],10) + current_from_flux(V[k]-inter_fermi,Eic[i],Eg[i]))
      P[k] = V[k] * J[k]
    max_eff_1[i] = np.amax(P)/power_in # Max efficiency for random Ei and Eg
  max_eff_2[j] = np.amax(max_eff_1) # Max efficiency for random Ei but correct Eg
  max_Eic[j] = Eic[np.argmax(max_eff_1)] # Record best Eic to report Eiv/Eic pairs
  print 'Best Eiv: {0} and Eic: {1}'.format(Eiv[j],max_Eic[j])
  
  
#Vmax[j] = optimize.fminbound(power_from_flux,0,E[j],args=(E[j],10))
  


plt.plot(Eiv, max_eff_2, 'o', c = 'b')
plt.show()