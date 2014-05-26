import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from scipy.optimize import newton

###
#		PARAMETERS FOR MOLECULE
###
molecule = 'methane'
Tc = -82.59 + 273.15 # Critical temperature, K
Pc = 45.99 # Critical pressure, bar
omega = 0.011 # accentric factor
R = 8.314e-5 # universal gas constant, m3-bar/K-mol
print "Molecule: %s. Critical Temp = %.1f K. Critical Pressure = %.1f bar." % (molecule, Tc, Pc)

###
#		GET CONDITIONS FROM INPUT
###
if len(sys.argv) != 3:
	sys.exit('Run program as:\npython PREOS.py T_in_K P_in_bar')
T = float(sys.argv[1])
P = float(sys.argv[2])

###
#		PREOS, BUILD CUBIC TO SOLVE FOR COMPRESSIBILITY
#   http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
###
Tr = T / Tc
a = 0.457235 * R**2 * Tc**2 / Pc
b = 0.0777961 * R * Tc /Pc
kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
alpha = ( 1 + kappa * (1 - np.sqrt(Tr) ) )**2

A = a * alpha * P / R**2 / T**2
B = b * P / R /T
## build cubic
def g(z):
	return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (A * B - B**2 - B**3)
z = newton(g,1) # compressibility factor
rho = P / (R * T * z) 

fugacity_coeff = np.exp(z - 1 - np.log(z - B) - A / np.sqrt(8) / B * np.log( (z + (1 + np.sqrt(2)) * B) / (z + (1 - np.sqrt(2)) * B) ))
print "Temperature %.2f K, Pressure %.2f bar" % (T, P)
print "Compressibility factor : %f" % z
print "Fugacity coefficient: %f" % fugacity_coeff
print "Fugacity at pressure %.3f bar = %.3f bar" % (P , fugacity_coeff * P)
print "Density: %f mol/m3" % ( rho )
print "Molar volume: %f L/mol" % ( 1/ rho * 1000 )
print "Density: %f v STP/v" % ( rho * 22.4 / 1000 )
print "IG density: %f v STP/v" % (rho * 22.4/ 1000 * z)
print "Energy density: %f MJ/L" % (rho *16.04 /1000.0/1000.0*55.5)

zz = np.linspace(0, 1.5)
fig = plt.figure()
plt.plot(zz , g(zz),color='k')
plt.xlabel('Compressibility')
plt.ylabel('g(z)')
plt.axvline(x = z )
plt.axhline(y = 0 )
plt.title('Visualizing roots')
plt.show()
