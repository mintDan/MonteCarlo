import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

a = 0
b = 10
n = 500
x = np.linspace(a,b,n)
#fplot = np.array(x**2)

def f(x):
	return x**2


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
print('Probability')
print(P)



Area = P*(b-a)*f(x[-1])
print('Approximated area')
print(Area)

print('Analytical area')
AArea = (1/3)*10**3
print(AArea)


plt.plot(x,f(x))
#pyplot.hold(True)
#pyplot.plot(x,f(x),c='black')
#pyplot.scatter(xg,yg,color='blue')
#pyplot.scatter(xg[inse],yg[inse],color='black')
#pyplot.xlabel('x')
#pyplot.ylabel('y')
#pyplot.title('Integrate $x^2$, area = {0:2f}'.format(AArea))
#pyplot.legend(['Approx Area = {0:2f}'.format(Area)])
#pyplot.savefig('MCint.png', bbox_inches='tight')
plt.show()