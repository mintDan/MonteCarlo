from matplotlib import pyplot
import numpy as np

#Grid
x0 = -2
x1 = 2
y0 = -2
y1 = 2
#Circle
a = -0.5
b = 0.5
r=np.sqrt(b**2)
n = 1000

#x = np.linspace(x0,x1,n)


#To plot the function
xr = np.linspace(a,b,n)
def f1(x):
	return np.sqrt(r**2-xr**2)
def f2(x):
	return -np.sqrt(r**2-xr**2)

xg = np.random.uniform(x0, x1, size=n)
yg = np.random.uniform(y0, y1, size=n)
# xg = np.random.uniform(a, b, size=n)
# yg = np.random.uniform(a, b, size=n)

inside = 0
for i in range(n):
	if xg[i]**2+yg[i]**2 <= r**2:
		inside+=1



		
P = inside/n
print('Probability')
print(P)


print('Analytical area')
print(np.pi*r**2)


# Pi = P*(b-a)**2/(r**2)
# print('Approximated Pi')
# print('Pi={}'.format(Pi))

Pi = P*(x1-x0)*(y1-y0)/(r**2)
print('Approximated Pi')
print('Pi={}'.format(Pi))

pyplot.hold(True)
pyplot.scatter(xg,yg,color='blue')
pyplot.plot(xr,f1(xr),color='black')
pyplot.plot(xr,f2(xr),color='black')
#pyplot.axis([x0,x1,y0,y1])
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.title('Approximate pi=3.14')
pyplot.savefig('MCpi.png', bbox_inches='tight')
pyplot.show()