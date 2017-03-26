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

def q(x):
	return (3.0/2.0)*np.sqrt(x)
	
def puniform(x):
	return 1
	
def qinverse(y):
	return (3.0/2.0)*1.0/(2.0*np.sqrt(y))

x = np.random.uniform(0,1,100000)
xdraw = q(x)
xdraw2 = qinverse(x)

f1 = f(x)
f2 = f(xdraw)
f3 = f(xdraw2)

N = len(f1)

xsorted = np.sort(x)	
dxs = np.zeros(N)

#Her kan gøres med numpy slicing
for i in range(1,N):
	dxs[i] = (xsorted[i]-xsorted[i-1])
dxs[0] = xsorted[0]	

print("Total points")
print(N)
print("Wiki way")



y = x
def xfromy(y):
	return y**(2.0/3.0)
xinv = xfromy(y)
f4 = f(xinv)



V = 1
integral1 = V*np.sum(f1)/N
print(integral1)

w = puniform(x)/q(x)
integral2 = V*np.sum(f2*w)/N
print(integral2)


w = q(xdraw2)
integral3 = V*np.sum(f3/w)/N
print(integral3)

w = q(xinv)
integral4 = V*np.sum(f4/w)/N
print(integral4)



print("Doing New integral")
def f5(x):
	#fs = np.zeros(len(x))
	#fs[x<=5] = x[x<=5]
	#fs[x<=0] = 0
	#fs[x>5] = -x[x>5]+10
	#fs[x>=10] = 0
	
	#return fs
	return np.array(x**2)



#############################
#Triangle distribution
def weightw(x):
	ws = np.zeros(len(x))
	ws[x<=5] = x[x<=5]
	ws[x>5] = -x[x>5]+10
	ws *= 1/25
	return ws

def xfromy(y):
	xs = np.zeros(len(y))
	#frac = 5.0/6.0
	frac = 0.5
	xs[y<=frac] = np.sqrt(50*y[y<=frac])
	xs[y>frac] = 10-25*np.sqrt((10/25)**2-(2/25)*(1+y[y>frac]))
	return xs
y = np.random.uniform(0,1,10000)
xinvtriangle = xfromy(y)
f5xinv = f5(xinvtriangle)
N = len(y)
wtriangle = weightw(xinvtriangle)
integralTriangle = np.sum(f5xinv/wtriangle)/N
print(integralTriangle)


########################
#Uniform distribution
def weightwuniform(x):
	return np.zeros(len(x))+1/10
	
xuni = np.random.uniform(0,10,10000)
N = len(xuni)
f5uni = f5(xuni)
V = 10
integralUniform = V*np.sum(f5uni)/N
print(integralUniform)


####################
#Gaussian distribution
def weightwgaussian(x):
	sigma = 2
	return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(x-5)**2/(2*sigma**2))

def xfromygaussian(y):
	xs = np.zeros(len(y))
	sigma = 2
	z = 2*y-1
	xs = 5+(np.sqrt(2)*sigma*np.sqrt(np.pi)/2)*(z+np.pi/12*z**3+(7*np.pi**2)/480*z**5+(127*np.pi**3)/40320*z**7
	+(4369*np.pi**4)/5806080*z**9+(34807*np.pi**5)/182476800*z**11)
	return xs

y = np.random.uniform(0,1,10000)
N = len(y)
xinvgauss = xfromygaussian(y)
f5gauss = f5(xinvgauss)

wgauss = weightwgaussian(xinvgauss)
integralGauss = np.sum(f5gauss/wgauss)/N

print(integralGauss)

###########################################
#Plotting
xplot = np.linspace(0,10,100)
wplot1 = 24*weightw(xplot)
wplot2 = 32*weightwuniform(xplot)
wplot3 = 28*weightwgaussian(xplot)
fplot = f5(xplot)



plt.figure()
ax = plt.gca()
#plt.plot(xplot,wplot1,color="red")
#plt.plot(xplot,wplot2,color="blue")
#plt.plot(xplot,wplot3,color="green")
#plt.plot(xplot,fplot,color="black")
plt.scatter(xinvtriangle[::200],24*wtriangle[::200],color="green")
plt.scatter(xuni[::200],32*weightwuniform(xuni)[::200],color="red")
plt.scatter(xinvgauss[::200],28*wgauss[::200],color="blue")
plt.legend(["Triangle {:0.1f}".format(integralTriangle),"Uniform {:0.1f}".format(integralUniform),"Gaussian {:0.1f}".format(integralGauss)])
plt.plot(xplot,fplot,color="black")
plt.axis([0,10,0,50])
ax.fill_between(xplot, fplot, np.zeros(len(fplot)),alpha=0.5,color="skyblue")
plt.xlabel("x")
plt.ylabel("y")


plt.title("MC integration, Importance Sampling")
#plt.scatter(xinvtriangle,f5xinv,color="black")
plt.savefig('MCis.png', bbox_inches='tight')
plt.show()

################################################################################
#https://www.wolframalpha.com/input/?i=int+sqrt(x)dx+from+0+to+1
#Ah vent, integral af sqrt(x)dx fra 0 til 1 er IKKE 1!!!!
#Det giver 2/3... ofc, fordi man får jo 2/3*x^3/2 når man integrate x^1/2... så yea...
# Der skal altså være 3/2 faktor på sqrt(x) distributioN!