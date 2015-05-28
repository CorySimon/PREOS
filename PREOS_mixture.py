import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from scipy.optimize import newton

###
#		PARAMETERS FOR MOLECULE
###
components = ['Xe','Kr']
Tc = np.array([16.59+273.15,-63.67+273.15]) # Critical temperature, K
Pc = np.array([58.42,55.25]) # Critical pressure, bar
# from cota
 # Tc = np.array([289.59,209.33])
 # Pc = np.array([5842000.0/100000.0,5525000.0/100000.0])
omega = np.array([0.0,0.0]) # accentric factor
R = 8.314e-5 # universal gas constant, m3-bar/K-mol
for i in range(2):
    print "Molecule: %s. Critical Temp = %.1f K. Critical Pressure = %.1f bar." % (components[i], Tc[i], Pc[i])

T = 298 # K 
P = 1.01325 # bar 
x = [0.2,0.8] # mole fraction
delta = - 0.0051 # binary interaction parameters (see reference in noblegas/readings)

###
#		PREOS, BUILD CUBIC TO SOLVE FOR COMPRESSIBILITY
#   http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
###
Tr = T / Tc
ac_i = 0.457235 * R**2 * Tc**2 / Pc
b_i = 0.0777961 * R * Tc /Pc
kappa_i = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
alpha_i = ( 1 + kappa_i * (1 - np.sqrt(Tr) ) )**2
a_i = alpha_i * ac_i

aij = (1-delta) * np.sqrt(a_i[0] * a_i[1])
# at this point, these are all 2D arrays. now use mixing rules
a = a_i[0] * x[0] ** 2 + a_i[1] * x[1] ** 2 + 2 * x[0] * x[1] * aij
b = x[0] * b_i[0] + x[1] * b_i[1]
# same as before
A = a * P / R**2 / T**2
B = b * P / R /T

## build cubic
def g(z):
	return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (A * B - B**2 - B**3)
Z = newton(g,1) # compressibility factor
rho = P / (R * T * Z) 

Lnfug_0 = -np.log(Z-B) + (Z-1)*b_i[0]/b - A / np.sqrt(8) / B *( 2.0/a * (x[0] *a_i[0] + x[1] * aij   ) -b_i[0]/b )*np.log( (Z+(1+np.sqrt(2))*B)/(Z+(1-np.sqrt(2))*B))
Lnfug_1 = -np.log(Z-B) + (Z-1)*b_i[1]/b - A / np.sqrt(8) / B *( 2.0/a * (x[0] *aij    + x[1] * a_i[1]) -b_i[1]/b )*np.log( (Z+(1+np.sqrt(2))*B)/(Z+(1-np.sqrt(2))*B))
print "Temperature %.2f K, Pressure %.2f bar" % (T, P)
print "Compressibility factor : %f" % Z
print "Fugacity coefficient, component %s: %f" % (components[0], np.exp(Lnfug_0))
print "Fugacity coefficient, component %s: %f" % (components[1], np.exp(Lnfug_1))
print "Fugacity component %s %.3f bar" % (components[0] ,P * np.exp(Lnfug_0) * x[0])
print "Fugacity component %s %.3f bar" % (components[1] ,P * np.exp(Lnfug_1) * x[1])
print "Density: %f mol/m3" % ( rho )
print "Molar volume: %f L/mol" % ( 1/ rho * 1000 )
print "Density: %f v STP/v" % ( rho * 22.4 / 1000 )
print "IG density: %f v STP/v" % (rho * 22.4/ 1000 * Z)
print "Energy density: %f MJ/L" % (rho *16.04 /1000.0/1000.0*55.5)

zz = np.linspace(0, 1.5)
fig = plt.figure()
plt.plot(zz , g(zz),color='k')
plt.xlabel('Compressibility')
plt.ylabel('g(z)')
plt.axvline(x = Z )
plt.axhline(y = 0 )
plt.title('Visualizing roots')
plt.show()
