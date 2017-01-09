"""
Using a simple monte carlo method to calculate a definite integral
modified 2017
Dan Krog 
"""
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt


def MC(a,b,n):
	"""
	This goes through the calculation of the definite integral
	Also calls the plotting function
	"""

	#Making array for plotting, doesn't have to be set with n.
	x = np.linspace(a,b,n)

	xg = np.random.uniform(a, b, size=n)
	yg = np.random.uniform(a, b**2, size=n)

	inse = np.array([],dtype=bool)
	below = 0
	for i in range(n):
		if yg[i] <= f(xg[i]):
			below+=1
			inse = np.append(inse,i)

	#int x^2 -> (1/3)*x^3
	#[(1/3)x^3]^1,0 = 1/3
	P = below/n
	print('Probability to land under function')
	print(P)



	Area = P*(b-a)*np.max(f(x))
	#print('Approximated area')
	#print(Area)

	#print('Analytical area')
	AArea = (1/3)*10**3
	#print(AArea)

	plot(x,f(x),xg,yg,xg[inse],yg[inse],AArea,Area)


def plot(x,fx,xg,yg,xginse,yginse,AArea,Area):
	plt.hold(True)
	plt.plot(x,fx,color='black')
	plt.scatter(xg,yg,color='blue')
	plt.scatter(xginse,yginse,color='black')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Integrate $x^2$, area = {0:2f}'.format(AArea))
	plt.legend(['Approx Area = {0:2f}'.format(Area)])
	plt.savefig('figs/MCint.png', bbox_inches='tight')
	plt.show()



if __name__ == "__main__":

	x0 = 0 #left boundary
	x1 = 10 #right boundary
	n = 500 #number of dots/random points

	#function to be integrated
	def f(x):
		return x**2

	MC(x0,x1,n)