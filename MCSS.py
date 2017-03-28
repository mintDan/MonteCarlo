"""
Using a simple monte carlo method to calculate a definite integral
modified 2017
Dan Krog 
"""
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata



def f(x):
	return np.array(x**2)

Nruns = 1000
Intabresults = np.zeros(Nruns)
Intcresults = np.zeros(Nruns)
	
N = 100000
Na = N/2
Nb = N/2

Va = 0.5
Vb = 0.5
Vc = 1
for i in range(Nruns):

	xa = np.random.uniform(0,0.5,int(Na))
	xb = np.random.uniform(0.5,1,int(Nb))
	xc = np.random.uniform(0,1,N)

	fa = f(xa)
	fb = f(xb)
	fc = f(xc)
	#fabavg = 0.5*(fa+fb)



	Intab = Va*np.sum(fa)/Na+Vb*np.sum(fb)/Nb
	Intabresults[i] = Intab
	#print("Integral of stratified sampling")
	#print(Intab)

	Intc = Vc*np.sum(fc)/N
	Intcresults[i] = Intc
	
	#print("Integral of crude Monte Carlo")
	#print(Intc)
	
	#Adaptive stratification
	avgfa = np.sum(fa)/Na
	sigmafa = np.sqrt((1/Na)*np.sum((fa-avgfa)**2))
	
	avgfb = np.sum(fb)/Nb
	sigmafb = np.sqrt((1/Nb)*np.sum((fb-avgfb)**2))
	
	Na = N*sigmafa/(sigmafb+sigmafa)
	Nb = N-Na
	

	
print("Average stratified")
Avgab = np.sum(Intabresults)/len(Intabresults)
print(Avgab)
print("Variance stratified")
print((1/len(Intabresults))*np.sum((Intabresults-Avgab)**2))

print("Average crude")
Avgc = np.sum(Intcresults)/len(Intcresults)
print(Avgc)
print("Variance Crude")
print((1/len(Intcresults))*np.sum((Intcresults-Avgc)**2))

plt.figure()
ax = plt.gca()
plt.scatter(xa,fa,c="black")
plt.scatter(xb,fb,c="black")
plt.axis([0,1,0,1])
#ax.fill_between(xplot, fplot, np.zeros(len(fplot)),alpha=0.5,color="skyblue")
plt.xlabel("x")
plt.ylabel("y")


plt.title("MC integration, Stratified Sampling")
#plt.scatter(xinvtriangle,f5xinv,color="black")

plt.show()


plt.figure(2)
ax2 = plt.gca()
plt.title("Histogram, integral of $x^2$")
#np.histogram(Intabresults)
#np.histogram(Intcresults)

plt.hist(Intcresults,alpha=0.5)
plt.hist(Intabresults,alpha=0.5)
plt.legend(["Crude","Stratified"])
plt.xlabel("Integral value")
plt.savefig('MCSS.png', bbox_inches='tight')
plt.show()