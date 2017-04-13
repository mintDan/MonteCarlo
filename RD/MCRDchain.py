import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.misc




Tend = 300
dt = 0.5
Nt = int(Tend/dt)
t = np.linspace(0,Tend,Nt+1)


NXstart = 1000
NYstart = 0
NXparticles = np.zeros(Nt+1)
NXparticles[0] = NXstart

NYparticles = np.zeros(Nt+1)
NYparticles[0] = NYstart

tauX = 7.2
tauY = 200


#Probability of decay, pr nucleus, 
#Hvis at Tend ikke er alt for stor, så kan man både skrive P = Lambda*Tend og P = Lambda*dt,
#Så skal se hvilken det er at sources mener specifically
PX = (1.0/tauX)*dt
PY = (1.0/tauY)*dt

#Halflife thallium i seconds
#T12 = 2007
#Lambda = np.log(2.0)/T12
#P = 1-2**(-dt/T12)



time = 0

for nt in range(1,Nt+1):

	NXparticlesnow = int(NXparticles[nt-1])

	DecayedParticlesX = 0
	if NXparticlesnow != 0:
	
		X = np.random.uniform(0,1,NXparticlesnow)
		#DecayedX = X[X<=P]
		DecayedParticlesX = len(X[X<=PX])
		
	NXparticles[nt] = NXparticles[nt-1]-DecayedParticlesX
	
	
	NYparticlesnow = int(NYparticles[nt-1])
	DecayedParticlesY = 0
	if NYparticlesnow != 0:
		Y = np.random.uniform(0,1,NYparticlesnow)
		
		DecayedParticlesY = len(Y[Y<=PY])
		
	NYparticles[nt] = NYparticles[nt-1] + DecayedParticlesX - DecayedParticlesY

	
	time+= dt
				
			
	

	


#=======================================
#Figure 1


plt.plot(t,NXparticles,color="black")
plt.plot(t,NYparticles,color="red")
plt.xlabel("time")
plt.ylabel("N particles")
plt.title("Monte Carlo simulation of Radioactive Decay chain, X->Y")
plt.legend(["X 210Bi","Y 210Po"])
plt.savefig('RDchain.png', bbox_inches='tight')
plt.show()

