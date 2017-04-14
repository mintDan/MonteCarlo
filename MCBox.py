import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.misc


Nparticles = 1000
nlstart = 1000


dt = 1
Tend = 4000
Ntime = int(Tend/dt)


t = np.linspace(0,Tend,Ntime+1)

nlarray = np.zeros(Ntime+1)
nlarray[0] = nlstart

#standard deviation stats
nlsquared = 0
nlavg = 0
count = 0

for nt in range(Ntime):
	
	ParticlesNow = int(nlarray[nt])
	
	Prob = ParticlesNow/Nparticles
	
	#X = np.random.uniform(0,1,ParticlesNow)
	#NToRight = len(X[X <= Prob])
	
	#Y = np.random.uniform(0,1,Nparticles-ParticlesNow )
	#NToLeft = len(Y[Y <= 1-Prob])
	Move = 0
	x = np.random.uniform(0,1)
	if x <= Prob:
		Move = Move - 1
	
	y = np.random.uniform(0,1)
	if y <= 1-Prob:
		Move = Move + 1
	
	nlarray[nt+1] = nlarray[nt] + Move
	
	if nt > Ntime/1.5:
		count += 1
		nlsquared += nlarray[nt+1]**2
		nlavg += nlarray[nt+1]
	#nlarray[nt+1] = nlarray[nt] - NToRight + NToLeft

nlsquared *=1/count
nlavg *=1/count
sigmanl = np.sqrt(nlsquared-nlavg**2)

print(sigmanl)


fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(121)

analyticalf = Nparticles/2*(1+np.exp(-2*t/Nparticles))

plt.plot(t,nlarray,color="black")
plt.plot(t,analyticalf,color="red")
plt.legend(["MC sim","Analytical sol"])
plt.title("Particles in a box converging to equilibrium state")
plt.xlabel("nt time")
plt.ylabel("N particles in left half")


x1 = np.random.uniform(0,0.45,150)
y1 = np.random.uniform(0,1,150)
x2 = np.random.uniform(0.55,1,30)
y2= np.random.uniform(0,1,30)

ax = fig.add_subplot(122)
plt.scatter(x1,y1,color="black")
plt.scatter(x2,y2,color="black")
plt.plot([0.5,0.5],[0,0.45],color="red",linewidth=5)
plt.plot([0.5,0.5],[0.55,1],color="red",linewidth=5)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Particle diffusion")

plt.savefig('Box.png', bbox_inches='tight')
plt.show()
