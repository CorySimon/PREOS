import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from scipy.optimize import newton

R = 8.314e-5  # universal gas constant, m3-bar/K-mol
class Molecule:
    """
    Store molecule info here
    """
    def __init__(self, name, Tc, Pc, omega):
        # name
        self.name = name
        # Critical temperature (K)
        self.Tc = Tc
        #! Critical pressure (bar)
        self.Pc = Pc
        #! Accentric factor
        self.omega = omega

    def print_params(self):
        print """Molecule: %s.
        \tCritical Temperature = %.1f K
        \tCritical Pressure = %.1f bar.
        \tAccentric factor = %f""" % (self.name, self.Tc, self.Pc, self.omega)

methane = Molecule("methane", -82.59 + 273.15, 45.99, 0.011)
methane.print_params()

def PREOS(molecule, T, P, plotcubic=False, printresults=True):
    """
    Peng-Robinson equation of state (PREOS)
    http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
    :param T float: temperature in Kelvin
    :param P float: pressure in bar
    """
    # build params in PREOS
    Tr = T / molecule.Tc  # reduced temperature
    a = 0.457235 * R**2 * molecule.Tc**2 / molecule.Pc
    b = 0.0777961 * R * molecule.Tc / molecule.Pc
    kappa = 0.37464 + 1.54226 * molecule.omega - 0.26992 * molecule.omega**2
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
    
    if printresults:
        print """PREOS calculation at
        \t T = %.2f K
        \t P = %.2f bar""" % (T, P)
        print "\tCompressibility factor : %f" % z
        print "\tFugacity coefficient: %f" % fugacity_coeff
        print "\tFugacity at pressure %.3f bar = %.3f bar" % (P , fugacity_coeff * P)
        print "\tDensity: %f mol/m3" % ( rho )
        print "\tMolar volume: %f L/mol" % ( 1/ rho * 1000 )
        print "\tDensity: %f v STP/v" % ( rho * 22.4 / 1000 )
        print "\tDensity of ideal gas at same conditions: %f v STP/v" % (rho * 22.4/ 1000 * z)
    
    if plotcubic:
        # Plot the cubic equation to visualize the roots
        zz = np.linspace(0, 1.5)  # array for plotting

        fig = plt.figure()
        plt.plot(zz , g(zz), color='k')
        plt.xlabel('Compressibility')
        plt.ylabel('g(z)')
        plt.axvline(x=z)
        plt.axhline(y=0)
        plt.title('Visualizing roots')
        plt.show()
    return rho # mol/m3

PREOS(methane, 298.0, 1.0)
