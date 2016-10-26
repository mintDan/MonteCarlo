# Monte Carlo simulations
Some scripts based on Monte Carlo methods.

## MC integration
![MCint.png](https://github.com/Bootlegg/MonteCarlo/blob/master/MCint.png)

## MC calculate pi
![MCpi.png](https://github.com/Bootlegg/MonteCarlo/blob/master/MCpi.png)

## Buffon's Needle
![MCbuffon.png](https://github.com/Bootlegg/MonteCarlo/blob/master/MCBuffon.png)

##Buffon's Needle Normal Distribution
Instead of using a uniform distribution to throw the needle, here are normal distributed needles. 
Pi can still be calculated explicitly, but an approximation to the probability is used in practice, since the normal distribution in essence goes to +- infinity.
The approximated probability P to land near the lines is given below

![P.png](https://github.com/Bootlegg/MonteCarlo/blob/master/P.png)

Where the integral I is given by  

![I.png](https://github.com/Bootlegg/MonteCarlo/blob/master/I.png)

The real value for P will be slightly larger than this, which the obtained values for pi show, since it always overshoots.

To calculate the probability, we need to calculate the integral I, and this in itself is done by changing the MonteCarlo.py script.

![BFint.png](https://github.com/Bootlegg/MonteCarlo/blob/master/BFint.png)

With these results we can continue to find pi, by isolating it from the expressions above.  

[pi.png](https://github.com/Bootlegg/MonteCarlo/blob/master/pi.png)

The result is shown below  

![MCbuffon.png](https://github.com/Bootlegg/MonteCarlo/blob/master/MCBuffonGauss.png)  
