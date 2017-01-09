"""Doing Buffon's Needle but with normal distribution of needles...
Calls the MonteCarlo module to calculate an integral
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import MonteCarlo



#Calculating needed integral by calling MonteCarlo module
x0 = 0
x1 = 0.1
n = 1000
def f(x):
	return (1.0/np.sqrt(2.0*0.2**2))*np.arccos(x/0.1)*np.exp(-(x-0.5)**2/(2.0*0.2**2))

integr = MonteCarlo.MC(x0,x1,n,f)



#Doing Buffon's Needle experiment with 
N = 2000000

x0 = 0
x1 = 1
y0 = 0
y1 = 1
l = 0.1
sigma = 0.2
mu = 0.5

#Throwing the needles

x = np.random.normal(mu, sigma, size=N)
y = np.random.uniform(y0, y1, size=N)

theta = np.random.uniform(0,np.pi, size=N)

#Checking if they're crossing the lines
n = 0
for i in range(N):
    if x[i]+l*np.cos(theta[i]) > x1:
        n+=1
    elif x[i]+l*np.cos(theta[i]) < x0:
        n+=1
#Experimental probability
P = n/N
print("Experimental probability")
print(P)

#Ptheo
#Ptheo2
#Without pi in the normalization factor
#integr = 0.025920
Ptheo = 4*integr/np.pi**(3.0/2.0)
print("Theoretical probability")
print(Ptheo)



pi = (4.0*integr/P)**(2.0/3.0)
print('pi={}'.format(pi))


#Let's make normal distribution function to plot
xp = np.linspace(x.min(),x.max(),100)
#with normalization (1/np.sqrt(2*np.pi*sigma**2))
nplot = np.array(0.33*np.exp(-(xp-mu)**2/(2*sigma**2)))+0.5

#Plotting 
#plt.hold(True)
fig = plt.figure()
ax = plt.gca()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Needles')
plt.axis([mu-5*sigma,mu+5*sigma,0-sigma,1+sigma])
for i in range(int(np.floor(N/10000))):
    ax.plot([x[i],x[i]+l*np.cos(theta[i])],[y[i],y[i]+l*np.sin(theta[i])],'b-')
ax.plot([x0,x0],[y0,y1],'-',color='black')
ax.plot([x1,x1],[y0,y1],'-',color='black')
ax.plot(xp,nplot,color='blue')
ax.legend(['$\pi$ = {0:2f}'.format(pi)])

fig.savefig('figs/MCBuffonGauss.png', bbox_inches='tight')
plt.show()



