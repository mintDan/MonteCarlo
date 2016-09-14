from matplotlib import pyplot
import numpy as np

a = 0
b = 10
n = 500
x = np.linspace(a,b,n)
#fplot = np.array(x**2)

def f(x):
	return x**2


xg = np.random.uniform(a, b, size=n)
yg = np.random.uniform(a, b**2, size=n)

below = 0
for i in range(n):
	if yg[i] <= f(x[i]):
		below+=1

#int x^2 -> (1/3)*x^3
#[(1/3)x^3]^1,0 = 1/3
P = below/n
print('Probability')
print(P)



Area = P*(b-a)*f(x[-1])
print('Approximated area')
print(Area)

print('Analytical area')
print((1/3)*10**3)


pyplot.hold(True)
pyplot.plot(x,f(x),c='black')
pyplot.scatter(xg,yg,color='blue')
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.title('Integrate $x^2$')
pyplot.savefig('MCint.png', bbox_inches='tight')
pyplot.show()