import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from scipy.optimize import newton


"""
Enter parameters for the molecule here
"""
molecule = 'methane'
Tc = -82.59 + 273.15  # Critical temperature, K
Pc = 45.99  # Critical pressure, bar
omega = 0.011  # accentric factor
R = 8.314e-5  # universal gas constant, m3-bar/K-mol
print "Molecule: %s.\n\tCritical Temperature = %.1f K\n\tCritical Pressure = %.1f bar." % (molecule, Tc, Pc)

"""
Take in arguments. Run in terminal as:
    python PREOS.py $temperature $pressure
   temperature in Kelvin 
   pressure in bar
"""
if len(sys.argv) != 3:
    sys.exit('Run program as:\npython PREOS.py T_in_K P_in_bar')
T = float(sys.argv[1])
P = float(sys.argv[2])

"""
Peng-Robinson equation of state: build cubic polynomial to solve for compressibility
http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
"""
# build params in PREOS
Tr = T / Tc  # reduced temperature
a = 0.457235 * R**2 * Tc**2 / Pc
b = 0.0777961 * R * Tc /Pc
kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2

A = a * alpha * P / R**2 / T**2
B = b * P / R /T

# build cubic polynomial
def g(z):
    return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (A * B - B**2 - B**3)

"""
Solve cubic polynomial for the compressibility factor
"""
z = newton(g,1)  # compressibility factor
rho = P / (R * T * z)  # density

# fugacity coefficient comes from an integration
fugacity_coeff = np.exp(z - 1 - np.log(z - B) - A / np.sqrt(8) / B * np.log( (z + (1 + np.sqrt(2)) * B) / (z + (1 - np.sqrt(2)) * B) ))

# print results
print "Temperature %.2f K, Pressure %.2f bar" % (T, P)
print "Compressibility factor : %f" % z
print "Fugacity coefficient: %f" % fugacity_coeff
print "Fugacity at pressure %.3f bar = %.3f bar" % (P , fugacity_coeff * P)
print "Density: %f mol/m3" % ( rho )
print "Molar volume: %f L/mol" % ( 1/ rho * 1000 )
print "Density: %f v STP/v" % ( rho * 22.4 / 1000 )
print "Density of ideal gas at same conditions: %f v STP/v" % (rho * 22.4/ 1000 * z)

"""
Plot the cubic equation to visualize the roots
"""
zz = np.linspace(0, 1.5)  # array for plotting

fig = plt.figure()
plt.plot(zz , g(zz), color='k')
plt.xlabel('Compressibility')
plt.ylabel('g(z)')
plt.axvline(x=z)
plt.axhline(y=0)
plt.title('Visualizing roots')
plt.show()
