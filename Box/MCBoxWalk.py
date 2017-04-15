import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.misc


Nparticle = 1000
dt = 1
#Tend = 40

x2 = 10
dx = 1
nx = int(x2/dx + 1)


D = dx**2/(4*dt)



Nparticlex = np.zeros(Nparticle)
Nparticley = np.array([y/Nparticle for y in range(Nparticle)])

#Nparticlex = np.random.randint(0,10,Nparticle)

#Nparticlesx = np.zeros((1,Nparticle))
#Nparticles[0,:] = np.random.randint(0,10,1000)
#print(Nparticlex)
#Skal plot dem op af y, basically..

#Er der mon egentlig basically 50% chance for left og right movement? Mon?
#Medmindre vi er langs edge?
#Eller laver man percentage på en anden måde?

#Hvis vi laver 1 timestep ad gangen fx, 1 timestep = 1 dx, så kan vi højest have 10 timesteps for sure

# time = 0
# timem = 0
# for m in range(Tend):
	# for i in range(Nparticle):
		# x = np.random.uniform(0,1)
		# walk = 0
		# #man kunne måske også stå still, så P = 1/3?
		# if x <= 0.25:
			# walk += 1
		# elif x >= 0.75:
			# walk += -1
			
		
		# if walk == 1:
			# if Nparticlex[i] != x2:
				# Nparticlex[i] += walk
				
		# if walk == -1:
			# if Nparticlex[i] != 0:
				# Nparticlex[i] += walk
	
	# time += dt
	# timem += 1
		

def MonteCarloDiffusion():
	for i in range(Nparticle):
		x = np.random.uniform(0,1)
		walk = 0
		#man kunne måske også stå still, så P = 1/3?
		if x <= 0.25:
			walk += 1
		elif x >= 0.75:
			walk += -1
			
		
		if walk == 1:
			if Nparticlex[i] != x2:
				Nparticlex[i] += walk
				
		if walk == -1:
			if Nparticlex[i] != 0:
				Nparticlex[i] += walk
	
	#time += dt
	#timem += 1
		

	Cx = np.zeros(nx)
	
	for i in range(Nparticle):
		Cx[int(Nparticlex[i])] += 1
		
	#print(sum(Cx))
	
	#xplot = np.linspace(0,Xend,Nx+1)
	#xnplot = np.array([n for n in range(11)])

	
	return Cx



#Figure object
fig = plt.figure()

#Subplot 1
ax1 = fig.add_subplot(121)
f = 10
ax1.scatter(Nparticlex[::f],Nparticley[::f],color="black")
ax1.set_xlabel("x")
ax1.set_ylabel("N")
ax1.set_title("Diffusion in 1D")
ax1.axis([0,10,0,1])


#Subplot 2
xnplot = np.linspace(0,x2,nx)
print(xnplot)
ax2 = fig.add_subplot(122)
ax2.legend(["MonteCarlo","Analytical"])
ax2.set_xlabel("x")
ax2.set_ylabel("concentration c")
ax2.axis([0,10,0,1000])

ax2.set_title("Concentration at given timestep")



def animate(i): #i increment with 1 each step
	ax1.clear()
	ax2.clear()
	
	#ax1.clf()
	#ax2.clf()
	
	#ax1.cla()
	#ax2.cla()
	#plt.close(fig)

	#plt.delaxes()
	#del ax1.lines
	#del ax2.lines
	#plt.cla()
	#del ax1.collections[:]
	#del ax2.collections[:]
	
	#try:
		#del plot1
		#del plot2
	#except:
		#pass

	
	if i > 0:
		Cxplot = MonteCarloDiffusion()
	
		#fig = plt.figure()

		#ax = fig.add_subplot(121)
		#f = 10
		ax1.scatter(Nparticlex[::f],Nparticley[::f],color="black")
		ax1.set_xlabel("x")
		ax1.set_ylabel("N")
		ax1.set_title("Diffusion in 1D")
		ax1.axis([0,10,0,1])


		#time = (Tend+0.05)*150
		#c0 = 3000
		#AnalyticalSol = c0/np.sqrt(4*np.pi*d*time)*np.exp(-xplot**2/(4*d*time))
		timem = i
		AnalyticalSol = 2*Nparticle/(dx*np.sqrt(np.pi*timem))*np.exp(-xnplot**2/timem)
		#print(i)
		#plt.figure(2)
		#ax = fig.add_subplot(122)
		plot1 = ax2.plot(xnplot,Cxplot,color="black")
		plot2 = ax2.plot(xnplot,AnalyticalSol,color="red")
		ax2.legend(["MonteCarlo","Analytical"])
		ax2.set_xlabel("x")
		ax2.set_ylabel("concentration c")
		ax2.axis([0,10,0,1000])
		ax2.set_title("Concentration at timestep {}".format(i))

		#plot1.remove()
		#plot2.remove()
		

	if i == 30:
		fig.savefig('BoxAnim.png', bbox_inches='tight')
	return None


anim = animation.FuncAnimation(fig, animate, interval=100)


plt.show()