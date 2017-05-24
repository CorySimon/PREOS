Peng-Robinson Equation of State solver
======================================

This Python module computes the fugacity coefficient, compressibility, and density of a gas at a particular temperature and pressure using the Peng-Robinson Equation of State.

The Peng-Robinson equation of state characterizes gas molecules by the following parameters:
* `Tc`: critical temperature (K)
* `Pc`: critical pressure (bar)
* `omega`: accentric factor

As an example calculation, we consider methane at 65.0 bar and 298.0 K. Methane has a critical temperature of -82.59 deg. C and a critical pressure of 45.99 bar. Its accentric factor is 0.011. We first create a methane molecule object and print its stored parameters:

```python
import preos
# pass name, Tc, Pc, omega
methane = preos.Molecule("methane", -82.59 + 273.15, 45.99, 0.011)
methane.print_params()
```

Finally, to print the properties of methane at 65.0 bar and 298.0 K using the Peng-Robinson equation of state:

```python
# pass the Molecule, T, P of interest
props = preos.preos(methane, 298.0, 65.0, plotcubic=True, printresults=True)
```

The function `preos` returns a dictionary of properties of the gas at this `T` and `P`. The `plotcubic` flag will plot the cubic polynomial in the compressibility factor for you to manually check that the correct solution was identified. By passing `printresults=False`, the printed output will be suppressed.

If one knows the fugacity and wants to know the corresponding pressure:
```python
# pass the Molecule, T, f of interest
props = preos.preos_reverse(methane, 298.0, 56.739, plotcubic=False, printresults=True)
```

The function `preos_reverse` returns a dictionary of properties of the gas at this `T` and `f`. The `plotcubic` flag will plot the cubic polynomial in the compressibility factor for you to manually check that the correct solution was identified. By passing `printresults=False`, the printed output will be suppressed.

For a binary mixture, we specify the temperature, total pressure `P_T`, and mole fractions `x`. We have an addition parameter, `delta`, the binary interaction parameter between the two gases. For example, for a Xe/Kr mixture:

```python
xe = preos.Molecule("Xe", 16.59 + 273.15, 58.42, 0.0)
kr = preos.Molecule("Kr", -63.67 + 273.15, 55.25, 0.0)

T = 298 # K 
P_total = 1.01325 # bar 
x = [0.2, 0.8] # mole fractions
delta = - 0.0051 # binary interaction parameter for Xe/Kr

props = preos.preos_mixture(xe, kr, delta, T, P_total, x)
```
