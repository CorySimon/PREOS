import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

R = 8.314e-5  # universal gas constant, m3-bar/K-mol
class Molecule:
    """
    Store molecule info here
    """
    def __init__(self, name, Tc, Pc, omega):
        """
        Pass parameters desribing molecules
        """
        #! name
        self.name = name
        #! Critical temperature (K)
        self.Tc = Tc
        #! Critical pressure (bar)
        self.Pc = Pc
        #! Accentric factor
        self.omega = omega

    def print_params(self):
        """
        Print molecule parameters.
        """
        print("""Molecule: %s.
        \tCritical Temperature = %.1f K
        \tCritical Pressure = %.1f bar.
        \tAccentric factor = %f""" % (self.name, self.Tc, self.Pc, self.omega))

def preos(molecule, T, P, plotcubic=True, printresults=True):
    """
    Peng-Robinson equation of state (PREOS)
    http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
    :param molecule: Molecule molecule of interest
    :param T: float temperature in Kelvin
    :param P: float pressure in bar
    :param plotcubic: bool plot cubic polynomial in compressibility factor
    :param printresults: bool print off properties

    Returns a Dict() of molecule properties at this T and P.
    """
    # build params in PREOS
    Tr = T / molecule.Tc  # reduced temperature
    a = 0.457235 * R**2 * molecule.Tc**2 / molecule.Pc
    b = 0.0777961 * R * molecule.Tc / molecule.Pc
    kappa = 0.37464 + 1.54226 * molecule.omega - 0.26992 * molecule.omega**2
    alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2

    A = a * alpha * P / R**2 / T**2
    B = b * P / R / T

    # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (
                A * B - B**2 - B**3)

    # Solve cubic polynomial for the compressibility factor
    z = newton(g, 1.0)  # compressibility factor
    rho = P / (R * T * z)  # density

    # fugacity coefficient comes from an integration
    fugacity_coeff = np.exp(z - 1 - np.log(z - B) - A / np.sqrt(8) / B * np.log(
                (z + (1 + np.sqrt(2)) * B) / (z + (1 - np.sqrt(2)) * B)))

    if printresults:
        print("""PREOS calculation at
        \t T = %.2f K
        \t P = %.2f bar""" % (T, P))
        print("\tCompressibility factor : ", z)
        print("\tFugacity coefficient: ", fugacity_coeff)
        print("\tFugacity at pressure %.3f bar = %.3f bar" % (
                P, fugacity_coeff * P))
        print("\tDensity: %f mol/m3" % rho)
        print("\tMolar volume: %f L/mol" % (1.0 / rho * 1000))
        print("\tDensity: %f v STP/v" % (rho * 22.4 / 1000))
        print("\tDensity of ideal gas at same conditions: %f v STP/v" % (
                rho * 22.4/ 1000 * z))

    if plotcubic:
        # Plot the cubic equation to visualize the roots
        zz = np.linspace(0, 1.5)  # array for plotting

        plt.figure()
        plt.plot(zz, g(zz), color='k')
        plt.xlabel('Compressibility, $z$')
        plt.ylabel('Cubic $g(z)$')
        plt.axvline(x=z)
        plt.axhline(y=0)
        plt.title('Root found @ z = %.2f' % z)
        plt.show()
    return {"density(mol/m3)": rho, "fugacity_coefficient": fugacity_coeff,
            "compressibility_factor": z, "fugacity(bar)": fugacity_coeff * P,
            "molar_volume(L/mol)": 1.0 / rho * 1000.0}

def preos_reverse(molecule, T, f, plotcubic=False, printresults=True):
    """
    Reverse Peng-Robinson equation of state (PREOS) to obtain pressure for a particular fugacity
    :param molecule: Molecule molecule of interest
    :param T: float temperature in Kelvin
    :param f: float fugacity in bar
    :param plotcubic: bool plot cubic polynomial in compressibility factor
    :param printresults: bool print off properties

    Returns a Dict() of molecule properties at this T and f.
    """
    # build function to minimize: difference between desired fugacity and that obtained from preos
    def g(P):
        """
        :param P: pressure
        """
        return (f - preos(molecule, T, P, plotcubic=False, printresults=False)["fugacity(bar)"])

    # Solve preos for the pressure
    P = newton(g, f)  # pressure

    # Obtain remaining parameters
    pars = preos(molecule, T, P, plotcubic=plotcubic, printresults=printresults)
    rho = pars["density(mol/m3)"]
    fugacity_coeff = pars["fugacity_coefficient"]
    z = pars["compressibility_factor"]

    return {"density(mol/m3)": rho, "fugacity_coefficient": fugacity_coeff,
            "compressibility_factor": z, "pressure(bar)": P,
            "molar_volume(L/mol)": 1.0 / rho * 1000.0}

# TODO: Implement mixture in object-oriented way as well
def preos_mixture(molecule_a, molecule_b, delta, T, P_total, x, plotcubic=True, printresults=True):
    """
    Peng-Robinson equation of state (PREOS) for a binary mixture
    http://en.wikipedia.org/wiki/Equation_of_state#Peng.E2.80.93Robinson_equation_of_state
    :param molecule_a: Molecule molecule 1 of interest
    :param molecule_b: Molecule molecule 2 of interest
    :param delta: binary interaction parameter between molecules a and b
    :param T: float temperature in Kelvin
    :param P_total: float total pressure in bar
    :param x: array mole fractions
    :param plotcubic: bool plot cubic polynomial in compressibility factor
    :param printresults: bool print off properties
    """
    # build arrays of properties
    Tc = np.array([molecule_a.Tc, molecule_b.Tc])
    Pc = np.array([molecule_a.Pc, molecule_b.Pc])
    omega = np.array([molecule_a.omega, molecule_b.omega])

    # build params in PREOS
    Tr = T / Tc  # reduced temperature
    a0 = 0.457235 * R**2 * Tc**2 / Pc
    b = 0.0777961 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    a = a0 * (1 + kappa * (1 - np.sqrt(Tr)))**2
    
    # apply mixing rules
    aij = (1.0 - delta) * np.sqrt(a[0] * a[1])
    a_mix = a[0] * x[0]**2 + a[1] * x[1]**2 + 2.0 * x[0] * x[1] * aij
    b_mix = x[0] * b[0] + x[1] * b[1]
    A = a_mix * P_total / R**2 / T**2
    B = b_mix * P_total / R / T

    # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (
                A * B - B**2 - B**3)

    # Solve cubic polynomial for the compressibility factor
    z = newton(g, 1.0)  # compressibility factor
    rho = P_total / (R * T * z)  # density
    
    Lnfug_0 = -np.log(z - B) + (z - 1.0) * b[0] / b_mix - A / np.sqrt(8) / B * (2.0 / a_mix * (x[0] * a[0] + x[1] * aij) - b[0] / b_mix) *\
        np.log((z + (1.0 + np.sqrt(2)) * B) / (z + (1.0 - np.sqrt(2)) * B))
    Lnfug_1 = -np.log(z - B) + (z - 1.0) * b[1] / b_mix - A / np.sqrt(8) / B * (2.0 / a_mix * (x[1] * a[1] + x[0] * aij) - b[1] / b_mix) *\
        np.log((z + (1.0 + np.sqrt(2)) * B) / (z + (1.0 - np.sqrt(2)) * B))

    # fugacity coefficient comes from an integration
    fugacity_coefs = np.exp(np.array([Lnfug_0, Lnfug_1]))

    if printresults:
        print("""PREOS calculation at
        \t T = %.2f K
        \t P, total = %.2f bar""" % (T, P_total))
        print("\tDensity: %f mol/m3" % rho)
        print("\tCompressibility factor : %f" % z)
        print("Component 0, %s:" % molecule_a.name)
        print("\tFugacity coefficient: %f" % fugacity_coefs[0])
        print("\tFugacity: %f bar" % (fugacity_coefs[0] * x[0] * P_total))
        print("Component 1, %s:" % molecule_b.name)
        print("\tFugacity coefficient: %f" % fugacity_coefs[1])
        print("\tFugacity: %f bar" % (fugacity_coefs[1] * x[1] * P_total))

    if plotcubic:
        # Plot the cubic equation to visualize the roots
        zz = np.linspace(0, 1.5)  # array for plotting

        plt.figure()
        plt.plot(zz, g(zz), color='k')
        plt.xlabel('Compressibility, $z$')
        plt.ylabel('Cubic $g(z)$')
        plt.axvline(x=z)
        plt.axhline(y=0)
        plt.title('Root found @ z = %.2f' % z)
        plt.show()
    return {"density(mol/m3)": rho, "fugacity_coefficients": fugacity_coefs,
            "compressibility_factor": z}
