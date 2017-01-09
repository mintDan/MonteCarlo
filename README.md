# Monte Carlo simulations
Some scripts based on Monte Carlo methods.

1. [MC integration](https://github.com/mintDan/MonteCarlo#mc-integration)
2. [Approximate Pi](https://github.com/mintDan/MonteCarlo#mc-calculate-pi)
3. [Buffon's Needle with normal distribution](https://github.com/mintDan/MonteCarlo#buffons-needle-normal-distribution)


## MC integration
![MCint.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCint.png)

## MC calculate pi
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

To calculate the probability P, we need to calculate the integral I, and this in itself is done by changing the MonteCarlo.py script.

![BFint.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/BFint.png)

With these results we can continue to find pi. The result is shown below

![MCbuffon.png](https://github.com/mintDan/MonteCarlo/blob/master/figs/MCBuffonGauss.png)
