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
![MCbuffon.png](https://github.com/Bootlegg/MonteCarlo/blob/master/MCBuffonGauss.png)
To calculate the probability, we need to calculate an integral, and this in itself is done by changing the MonteCarlo.py script to calculate said integral.
![BFint.png](https://github.com/Bootlegg/MonteCarlo/blob/master/BFint.png)