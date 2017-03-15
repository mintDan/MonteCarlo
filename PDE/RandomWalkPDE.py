import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d


nx = 15
ny = 15
x = np.linspace(0,1,nx)
y = np.linspace(0,1,ny)

X,Y = np.meshgrid(x,y)


#I will basically do a random walk, ye?
#So, with probability something, we go either up,down,left,right...
#I need boundary conditions though

#En second order PDE i 2 dimensions, har 2x2 = 4 initial boundary conditions maybe?
#Så lad os sige

U = np.zeros((ny,nx))
U[0,:] = 0
U[:,0] = 0
U[ny-1,:] = 1 #Top row, right?
U[:,nx-1] = 0


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

nruns = 100000
#nsteps = 2000
for n in range(nruns):
	
	#Here we pick a point on the grid
	i = np.random.randint(0,nx-1)
	j = np.random.randint(0,ny-1)
	

	
	#point = [j,i]
	#we're at point U[j,i]
	
	
	#I can either do a frequency/count variable for each point, but atm the same point will just
	#be appended multiple times to the list, so in the end it works out as frequency anyway
	VisitedPoints = []
	VisitedUs = []
	#Here, we basically do the random walking of the point, so basically a loop in a loop
	
	#This is the current point we're at in the random walk
	jnow = j
	inow = i
	
	#This should be until we hit a boundary pretty much... but i think 10000 steps is good enough
	#for step in range(nsteps):
	OnBoundary = False
	while OnBoundary == False:
	
		if BoundaryTest(jnow,inow) == True:
			
			VisitedPoints.append([jnow,inow])
			VisitedUs.append(U[jnow,inow])
			
			#we break this inner loop if we arrived at boundary
			#break	
			OnBoundary = True
		
		else:
			P = np.random.uniform(0,1)
			if P <= 0.25:
				#We go to U[j,i+1]
				VisitedPoints.append([jnow,inow+1])
				VisitedUs.append(U[jnow,inow+1])
				inow += 1
			elif P > 0.25 and P <= 0.5:
				#we go to U[j+1,i]
				VisitedPoints.append([jnow+1,inow])
				VisitedUs.append(U[jnow+1,inow])
				jnow += 1
			elif P > 0.5 and P <= 0.75:
				#we go to U[j,i-1]
				VisitedPoints.append([jnow,inow-1])
				VisitedUs.append(U[jnow,inow-1])
				inow += -1
			else:
				#we go to U[j-1,i]
				VisitedPoints.append([j-1,i])
				VisitedUs.append(U[j-1,i])
				jnow += -1
	
	
		
	#So, we now start from the new point...
	
	#
	#Here, we do the averaging of the visited points to get U[j,i]
	#Oh wait... atm...
	#We RESET Us... we need to consider what if they weren't 0 already!
	#Right now, i take the average of what it was before, and what it is now..
	#That's probably not really correct either, though
	
	#Well, you can do it another way though...
	#The Tour Du Wino only records boundary points, apparently..
	#And in that case, we don't have to continually update
	#I think tuor du wino basically just for each point in the grid, it does n amount of tours...
	#But it doesn't pick the grid points randomly...
	#Right now, i pick randomly, and therefor, i might pick the same one twice, overwriting the previous result
	#if we systemically, with 2 nested loops, go through each point on a journey, then we don't overwrite...
	U[j,i] = (U[j,i] + np.sum(VisitedUs)/len(VisitedUs))/2
	
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Random Walk For Solving PDE")
ax.plot_surface(X, Y, U,cmap = 'jet')#,rstride=10, cstride=10)

#https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
#Jakevdp siger at jet er standard colormap, fra MatLab? Men den er bad?

plt.show()