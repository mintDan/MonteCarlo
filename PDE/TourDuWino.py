import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d

nruns = 100
nx = 25
ny = 25
Lx = 1
Ly = 1
dx = Lx/(nx-1)
dy = Ly/(ny-1)
x = np.linspace(0,Lx,nx)
y = np.linspace(0,Ly,ny)

X,Y = np.meshgrid(x,y)

U = np.zeros((ny,nx))
U[0,:] = 0
U[:,0] = 0
U[ny-1,:] = 1 #Top row
U[:,nx-1] = 0



#Let's look at uxx + sin(x)uyy = 0

#P1 = (1.0/4.0)*(1.0/(1+dx*dx*np.sin(X[j,i])/dy*dy))

def BoundaryTest(j,i):
	
	IsBoundary = False
	if j == 0:
		IsBoundary = True
	elif j == ny-1:
		IsBoundary = True
	elif i == 0:
		IsBoundary = True
	elif i == nx-1:
		IsBoundary = True
		
	return IsBoundary

for j in range(ny):
	for i in range(nx):
		VisitedUs = []
		
		for n in range(nruns):
		
			
			#Now we will systematically go through each point here, and do maybe 100 random walks for each
			#VisitedPoints = []
			
			#Here, we basically do the random walking of the point, so basically a loop in a loop
			
			#This is the current point we're at in the random walk
			jnow = j
			inow = i
			
			#This should be until we hit a boundary pretty much... but i think 10000 steps is good enough
			#for step in range(nsteps):
			OnBoundary = False
			while OnBoundary == False:
			
				if BoundaryTest(jnow,inow) == True:
					
					#VisitedPoints.append([jnow,inow])
					VisitedUs.append(U[jnow,inow])
					
					#we break this inner loop if we arrived at boundary
					#break	
					OnBoundary = True
				
				else:
					P = np.random.uniform(0,1)
					
					P1 = (1.0/2.0)*np.sin(X[j,i])/(dy*dy/(dx*dx)+np.sin(X[j,i]))
					
					P2 = (1.0/2.0)*(1.0/(1+dx*dx*np.sin(X[j,i])/(dy*dy)))
					P3 = P1
					
					if P <= P1:#0.25:
					
					
						#We go to U[j,i+1]
						#VisitedPoints.append([jnow,inow+1])
						#VisitedUs.append(U[jnow,inow+1])
						inow += 1
					elif P > P1 and P <= P1+P2:
						#we go to U[j+1,i]
						#VisitedPoints.append([jnow+1,inow])
						#VisitedUs.append(U[jnow+1,inow])
						jnow += 1
					elif P > P1+P2 and P <= P1+P2+P3:
						#we go to U[j,i-1]
						#VisitedPoints.append([jnow,inow-1])
						#VisitedUs.append(U[jnow,inow-1])
						inow += -1
					else:
						#we go to U[j-1,i]
						#VisitedPoints.append([j-1,i])
						#VisitedUs.append(U[j-1,i])
						jnow += -1
		
		U[j,i] = np.sum(VisitedUs)/len(VisitedUs)

		
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("U")
ax.set_title("Random Walk For PDE")
ax.plot_surface(X, Y, U,cmap = 'winter')#,rstride=10, cstride=10)

#https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
#Jakevdp siger at jet er standard colormap, fra MatLab? Men den er bad?
fig.savefig('TourDuWino.png', bbox_inches='tight')
plt.show()