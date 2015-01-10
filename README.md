Peng-Robinson Equation of State solver
======================================

This script computes the fugacity coefficient, compressibility, and density of a gas at a particular temperature and pressure using the Peng-Robinson Equation of State.

For your molecule, modify the script by entering the relevant parameters:
* `Tc`: critical temperature (K)
* `Pc`: critical pressure (bar)
* `omega`: accentric factor

Then, as an example, to find the properties at 298.0 K and 65.0 bar, run the script in the terminal as:
`python PREOS.py 298.0 65.0`

At the end, a plot will show the cubic polynomial to visualize the roots as a check that the correct solution was identified.
