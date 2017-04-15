# Monte Carlo simulations
Some scripts based on Monte Carlo methods.

1. [MC integration](https://github.com/mintDan/MonteCarlo#mc-integration)
2. [Ising model](https://github.com/mintDan/MonteCarlo#ising-model)
3. [Random Walk method](https://github.com/mintDan/MonteCarlo#random-walk-method)
4. [Radioactive Decay](https://github.com/mintDan/MonteCarlo#radioactive-decay)
5. [Particles in a box](https://github.com/mintDan/MonteCarlo#particles-in-a-box)
6. [Approximate Pi](https://github.com/mintDan/MonteCarlo#mc-calculate-pi)
7. [Buffon's Needle with normal distribution](https://github.com/mintDan/MonteCarlo#buffons-needle-normal-distribution)


## MC integration
The MonteCarlo.py script is a module that can be called from other scripts, and is used in Buffon's needle with normal distributed needles.
The probability is calculated by comparing areas here, so in theory only linearly distributed points can be used, I think. To use different distributions, one must 
use the traditional method as adapted from standard numerical integration.  

The script MCint3.py uses importance sampling with distributions seen in the legend, and the result of their integration. For this function, the integral is 333.33... analytically.


![MCint.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCis.png)

Another way of variance reduction is by stratified sampling. Even if you only split twice and give half the amount of points to each stratum, you'll still learn more about the integral compared to
during a crude MC integration. The cost of gaining the extra knowledge, despite keeping same points N in total, and giving N/2 to each stratum, is some extra computation. No free meal, but pretty close!

![MCSS.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCSS.png)

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


## Random Walk method
A partial differential equation can be solved by Random Walk method. Starting from some point, we do N random walks, where each walk ends when we land on a boundary. The value at the boundary is then recorded.
The frequency of boundaries can then give an estimate for the point, i.e if a point is close to a certain boundary, we expect a random walk from that point will end up at that boundary more often than a boundary further away.
So in a sense, the frequencies give weighted averages to the value of the point.

![TourDuWino.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/TourDuWino.png)

### Random Walk On Spheres
Instead of doing tiny steps while walking, we can make a dynamic stepping size by drawing a circle with size of radius as the distance to the nearest boundary point. The random walk will then be by going to a point on this drawn circle,
and checking if the new point is close to a boundary point. That way we take much larger steps, without sacrificing the essence of the Random Walk.
 We can never really expect the random walk to land on a boundary point, since our stepping is continous. Therefore we say that when the walk goes close enough, as dictated by a small epsilon,
the walk ends.
The added bonus is that this method works without grids and with arbitrary boundaries, and thus makes it easier to solve PDEs on arbitrary domains.
The epsilon that is used for detection when the walk is close to a boundary should probably scale with the distance between boundary points... And also, it happens that the walk might actually step OUTSIDE the domain, since the boundary in a sense is diffuse.

![WalkOnSpheres.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/WalkOnSpheresSubplot.png)

## Radioactive Decay
Simulating radioactive decay with MC method. The decay constant with the time interval gives the probability that a nucleus will decay, but it is really an approximation for small enough time intervals.  
With a short enough time step and large amount of initial nuclei, the curve will fitter to the theoretical curve better, and once we reach low amount of cores, the stochastic nature will show itself, since the exponential solution is really only an approximation.
For the first timestep, the expected number of decays is constant, so taking the number of decays from 1000 experiments, we see a poisson distributions of the decays, since the probability of decay is small.

![RD.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCdecay.png)

We can also simulate a chain reaction, nuclei X -> nuclei Y, where nuclei Y also are unstable. Starting with X = 1000 and Y = 0, we have a system of coupled differential equations.

![RDchaun.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/RDchain.png)


## Particles in a box
We can simulate diffusion by statistical methods, instead of keeping track of a big phase space, the positions and velocities of each particle, and calculating collisions, momentum conservations etc. Starting with all particles in the left half, there is a probability that the particles will go to the right.
Likewise, for particles in the right half, they have a probability to go to the left half. The probabilities will converge to each other, and we reach an equilibrium state.

![Box.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/Box.png)

### Animation for 1D diffusion, comparison MonteCarlo with analytical solution

The script MCboxwalk.py simulates diffusion of particles with a random walk. Each point has a 25% chance to go either left or right, and 50% chance to stay in place. The end result is much like what is seen from particle diffusion.

![BoxAnim.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/BoxAnim.png)

## MC calculate pi
Compares the amount of points landing inside and outside the circle, which together with the area outside, gives an approximation for pi.
As seen here the result is pi = 2.9.

![MCpi.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCpi.png)

## Buffon's Needle Normal Distribution
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
