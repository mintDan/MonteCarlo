import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.misc




Tend = 100
dt = 1
Nt = int(Tend/dt)
t = np.linspace(0,Tend,Nt+1)




Nstart = 200
Nparticles = np.zeros(Nt)
Nparticles[0] = Nstart

Lambda = 0.01


#Probability of decay, pr nucleus, 
#Hvis at Tend ikke er alt for stor, så kan man både skrive P = Lambda*Tend og P = Lambda*dt,
#Så skal se hvilken det er at sources mener specifically
P = Lambda*dt

#Halflife thallium i seconds
#T12 = 2007
#Lambda = np.log(2.0)/T12
#P = 1-2**(-dt/T12)

print("We must have BT << 1 for a good Poisson distribution")
print(Lambda*Nstart*dt)

print("Probability taylor approx must have lambda*t_interval << 1")
print(P)


Nexperiments = 1000

#DecayValues = np.zeros(Nt*Nexperiments)
DecayValues = np.zeros((Nt)*Nexperiments,dtype=int)
DecayValuesSum = np.zeros(Nexperiments,dtype=int)

#FOr a single experiment
#DecayValuesSingle = np.zeros(Nt)




#Measured Half Lives
HalfLives = np.zeros(Nexperiments)
HalfLivesSquared = np.zeros(Nexperiments)

for nexp in range(Nexperiments):
	
	Nparticles = np.zeros(Nt+1)
	Nparticles[0] = Nstart
	
	time=0
	
	
	#Testing for Thalf
	ReachedHalf = False
	#Thalf = 0
	
	
	DecayValuesSingle = np.zeros(Nt)
	
	for nt in range(1,Nt+1):

		Nparticlesnow = int(Nparticles[nt-1])
		
		DecayedParticles = 0
		if int(Nparticlesnow) != 0:
			#for npart in range(Nparticlesnow):
			#	X = np.random.uniform(0,1)
				
				#Seems like it should be less than a, if a is really small...
				#Or maybe not, dunno
				#It should be less than a, since a is probability of decay...
			#	if X =< P:
			#		DecayedParticles += 1
			
			X = np.random.uniform(0,1,Nparticlesnow)
			#DecayedX = X[X<=P]
			DecayedParticles = len(X[X<=P])
			
		
			#print(DecayedParticles)
			
		Nparticles[nt] = Nparticles[nt-1]-DecayedParticles
		
		#nt-1 fordi loop starter fra 1?
		if nt == 1:
			#print(nt)
			DecayValues[(nexp+1)*(nt)-1] = int(DecayedParticles)
			#print(DecayedParticles)
			DecayValuesSingle[nt-1] = int(DecayedParticles)
		#Single experiment, nt-1 fordi at for loop starter fra 1
		#DecayValuesSingle[nt-1] = int(DecayedParticles)
		
		time+= dt
		
		#Test for T1/2
		if ReachedHalf == False:
			if Nparticles[nt] <= Nstart/2:
				ReachedHalf = True
				#Thalf = time
				HalfLives[nexp] = time
				HalfLivesSquared[nexp] = time*time
				
	DecayValuesSum[nexp] = np.sum(DecayValuesSingle)	
	
print("Expected number of decays in a dt")
print(Nstart*dt*Lambda)

print("Maximum found decay in a dt")
print(max(DecayValues))

print("Measured avg Half life T1/2")
Thalfavg = np.average(HalfLives)
print(Thalfavg)
print("Theoretical Halftime T1/2")
print(np.log(2)/Lambda)
	
print("Sigma of half life")
AvgHalfLives = np.average(HalfLives)
AvgHalfLivesSquared = np.average(HalfLivesSquared)
sigmaThalf = np.sqrt((1/(Nexperiments-1))*(AvgHalfLives**2+AvgHalfLivesSquared))
print(sigmaThalf)
	
	
#=======================================
#I have the decay value arrays, maybe I can make a histogram out of them here, and not in the for loop

BinScatterY = np.bincount(DecayValuesSum)
BinScatterX = [int(x) for x in range(len(BinScatterY))]



#=======================================
#Figure 1
Ntheo = Nstart*np.exp(-Lambda*t)
fig = plt.figure(figsize=plt.figaspect(0.5))

#fig.suptitle("MC", fontsize=16)

ax = fig.add_subplot(121)

plt.plot(t,Nparticles,color="black")
plt.plot(t,Ntheo,color="red")
plt.xlabel("time")
plt.ylabel("N particles")
plt.title("Monte Carlo simulation of Radioactive Material")
plt.legend(["Monte Carlo","Theoretical"])
#plt.savefig('RD.png', bbox_inches='tight')
#plt.show()



#==========================================
#Figure 4, scatter

#Nx = [int(x) for x in range(int(DecayValuesSum.min()),int(DecayValuesSum.max()))]
Nx = np.linspace(int(DecayValuesSum.min()),int(DecayValuesSum.max()),20)
mu = Lambda*Nstart*dt
#Dtheo = Nstart*np.array([(mu**nx)*np.exp(-mu)/scipy.misc.factorial(nx) for nx in Nx])
Dtheo = Nexperiments*(mu**Nx)*np.exp(-mu)/scipy.misc.factorial(Nx)

ax = fig.add_subplot(122)
#plt.scatter(BinScatterX,BinScatterY,color="black")

plt.bar(BinScatterX,BinScatterY,color="black")
plt.plot(Nx,Dtheo,color="red")
plt.title("Poisson distribution of Decays")
plt.xlabel("Decays")
plt.ylabel("Counts")
plt.legend(["Expected Poisson","Monte Carlo test"])
plt.savefig('MCdecay.png', bbox_inches='tight')

plt.show()



#=======================================
#Figure 2

#print("Let's look at DecayValues")
#print(DecayValues)

Nx = [int(x) for x in range(int(DecayValues.min()),int(DecayValues.max()))]
#mu = Lambda*Nstart*dt
#mu = 4
#Dtheo = Nstart*np.array([(mu**nx)*np.exp(-mu)/scipy.misc.factorial(nx) for nx in Nx])

#DecaysTheoretical for large N, small p
sigma = np.std(DecayValues,ddof=1)
avg = np.average(DecayValues)
Dtheo = (1.0/np.sqrt(2*sigma**2*np.pi))*np.exp(-(Nx-avg)**2/(2*sigma**2))

#Kan måske også tage Nstart som scale?
plt.figure(3)
#plt.hist(DecayValues,normed=True)
plt.hist(DecayValues,normed=True)
#plt.hist(DecayValuesSingle)
plt.plot(Nx,Dtheo,color="red")
plt.xlabel("Ndecays")
plt.ylabel("Count")
plt.title("Poisson distribution of decays")
plt.savefig('Decays.png', bbox_inches='tight')
plt.show()


#==========================================
#Figure 3

plt.figure(4)
plt.hist(DecayValuesSum)#,normed=True)
plt.show()

