Peng-Robinson Equation of State solver
======================================

This Python module computes the fugacity coefficient, compressibility, and density of a gas at a particular temperature and pressure using the Peng-Robinson Equation of State.

The Peng-Robinson equation of state characterizes gas molecules by the following parameters:
* `Tc`: critical temperature (K)
* `Pc`: critical pressure (bar)
* `omega`: accentric factor

As an example calculation, we consider methane at 65.0 bar and 298.0 K. Methane has a critical temperature of -82.59 deg. C and a critical pressure of 45.99 bar. Its accentric factor is 0.011. We first create a methane molecule object and print its stored parameters:

    import preos
    # pass name, Tc, Pc, omega
    methane = preos.Molecule("methane", -82.59 + 273.15, 45.99, 0.011)
    methane.print_params()

Finally, to print the properties of methane at 65.0 bar and 298.0 K using the Peng-Robinson equation of state:
    
    # pass the Molecule, T, P of interest
    rho = preos.preos(methane, 298.0, 65.0, plotcubic=True, printresults=True)
    # returns density in mol/m3

The `plotcubic` flag will plot the cubic polynomial in the compressibility factor for you to manually check that the correct solution was identified. By passing `printresults=False`, the printed output will be supressed.
