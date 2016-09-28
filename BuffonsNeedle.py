import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



N = 1000000

x0 = 0
x1 = 1
y0 = 0
y1 = 1
l = 0.1

#Throwing the needles
x = np.random.uniform(x0, x1, size=N)
y = np.random.uniform(x0, x1, size=N)
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
print(P)

#Calculate pi
pi = 2*l/(P*(x1-x0))
print('pi={}'.format(pi))


#Plotting 
fig = plt.figure()
ax = plt.gca()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Needles')
for i in range(int(np.floor(N/10000))):
    ax.plot([x[i],x[i]+l*np.cos(theta[i])],[y[i],y[i]+l*np.sin(theta[i])],'b-')
ax.plot([x0,x0],[y0,y1],'-',color='black')
ax.plot([x1,x1],[y0,y1],'-',color='black')
ax.legend(['$\pi$ = {0:2f}'.format(pi)])
fig.savefig('MCBuffon.png', bbox_inches='tight')
plt.show() #fig.show() does not much?



