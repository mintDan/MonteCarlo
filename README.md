# Monte Carlo simulations
Some scripts based on Monte Carlo methods.

1. [MC integration](https://github.com/mintDan/MonteCarlo#mc-integration)
2. [Ising model](https://github.com/mintDan/MonteCarlo#ising-model)
3. [Approximate Pi](https://github.com/mintDan/MonteCarlo#mc-calculate-pi)
4. [Buffon's Needle with normal distribution](https://github.com/mintDan/MonteCarlo#buffons-needle-normal-distribution)


## MC integration
The MonteCarlo.py script is a module that can be called from other scripts, and is used in Buffon's needle with normal distributed needles.
The probability is calculated by comparing areas here, so in theory only linearly distributed points can be used, I think. To use different distributions, one must 
use the traditional method as adapted from standard numerical integration. 

![MCint.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCint.png)

## Ising model
This 2D Ising model sets up a square lattice with each site representing an electron spin with direction up +1 or down -1. Ferromagnetic behavior shows when there's a  
majority of spin either up or down, which can happen at a low temperature. At high temperature the thermal disturbances will randomize the direction of the magnetic dipoles(electron spin)  
such that the magnet becomes a paramagnet, losing its ferromagnetic properties.

![Ising.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Ising2.45.png)  

We see here the spin lattice at T = 2.45, slightly higher than the critical temperature. Domains are still visible, with red being spin up, blue being spin down.
Started from random configuration.

![Ising.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Ising4.85.png)  

For T = 4.85, the thermal disturbances are too strong for any significant domains to develop.


![Magnetization.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Magnetization.png)  

The magnetization is a measure of alignment and order, and therefore strength of the magnetic field from the magnet. Magnetization can be negative or positive depending on which way the macroscopic field is oriented.
Magnization is higher at low T, if you start from an ordered state, than at high T, since the thermal disturbances will break alignment.


![Energy.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Energy.png)  

Energy pr spin site. At low T, starting from an ordered state, the dipoles will try to be put into the lower energy states, as specified by the energy from their neighbors.


![Susceptibility.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Susceptibility.png)  

The susceptibility is a measure for how the material reacts to an external magnetic field. If x_v is positive, the material is more paramagnetic, and the dipoles inside can align to a greater degree
with the external field. The models gives low susceptibility for low T, since the dipoles are firm in their arrangement. Higher T will make the dipoles more free to be influenced by an external field, but,
too high T means the thermal disturbances will be too strong to give any significant alignment with an external field.


![SpecificHeat.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Specificheat.png)  

The specific heat describes how much the internal energy changes when the temperature changes. We see that the internal energy changes mostly around the critical temperature, which can also be seen from the energy graph.
Once the magnetization reaches 0, and the energy reaches it's maximum, we basically have random configuration of spins, and increasing T cannot do anything to increase the Energy in this system.


## MC calculate pi
Compares the amount of points landing inside and outside the circle, which together with the area outside, gives an approximation for pi.
As seen here the result is pi = 2.9.

![MCpi.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCpi.png)

##Buffon's Needle Normal Distribution
Instead of using a uniform distribution to throw the needle, here are normal distributed needles. 
Pi can still be calculated explicitly, but an approximation to the probability is used in practice, since the normal distribution in essence goes to +- infinity.
The approximated probability P to land near the lines is given below

![P.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/P.png)

The real value for P will be slightly larger than this, which the obtained values for pi show, since it always overshoots.

![Pi.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Pi.png)

Where the integral I is given by

![I.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/I.png)

To calculate the probability P, we need to calculate the integral I, and this in itself is done by calling MonteCarlo.py as a module.

![BFint.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/BFint.png)

With these results we can continue to find pi. The result is shown below

![MCbuffon.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCBuffonGauss.png)
