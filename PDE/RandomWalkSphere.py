import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d
from matplotlib.mlab import griddata

nx = 15
ny = 15
Lx = 1
Ly = 1
x = np.linspace(0,Lx,nx)
y = np.linspace(0,Ly,ny)

#Tolerance for hitting the edge
eps = 0.05



Boundary = []

for i in range(nx):
	#We do x,y,U, with 1 at y=1 and 0 everywhere else
	Boundary.append([x[i],0,0])
	Boundary.append([x[i],1,1])
	Boundary.append([0,y[i],0])
	Boundary.append([1,y[i],0])

#print(Boundary)
	
npoints = 75
nwalks = 100
nsteps = 2000

# will have data [x,y,U(x,y)]
Upoints = []

for n in range(npoints):
	
	#Here we pick a point on the grid
	xrandom = np.random.uniform(0,Lx)
	yrandom = np.random.uniform(0,Ly)
	
	
	FoundBoundaries = []#np.array([])
	#Either loops or nested function maybe
	

	for nwalk in range(nwalks):
		
		xnow = xrandom
		ynow = yrandom
		
		FoundBoundary = False
		while FoundBoundary == False:
		
			#Måske prøv at plotte dem... like, append disse points til en list, og plot dem
			#For at se om de kommer MEGET FAR AWAY from boundary, eller hvad der lige foregår...
			#Måske er det noget floating point error
			#if xnow > 1 or xnow < 0 or ynow > 1 or ynow < 0:
			#	print("Point is outside boundary!")
		
		
		#We check first if current point is near/at boundary
			#FoundBoundary = False
			for Bpoint in Boundary:
				d = np.sqrt((xnow-Bpoint[0])**2+(ynow-Bpoint[1])**2)
				
				if d < eps:
					#We found a boundary, and we can end this for loop and later the walk
					FoundBoundaries.append(Bpoint)
					FoundBoundary = True

					#We're also gonna break the loop, there's no need to go through ALL the rest of the points
					break
			
			if FoundBoundary == False:
				#Find biggest allowable radius, is also the minimum distance to boundary
				dmin = Lx
				for Bpoint in Boundary:
					d = np.sqrt((xnow-Bpoint[0])**2+(ynow-Bpoint[1])**2)
					if d < dmin:
						dmin = d
				
				#We go to a new point
				anglerandom = np.random.uniform(0,2*np.pi)
				
				#print(anglerandom)
				#New point is xnew = xold + dx, ynew = yold + dy
				xnow = xnow + dmin*np.cos(anglerandom)
				ynow = ynow + dmin*np.sin(anglerandom)
				

				
			#else:
				
				#reset xnow,ynow to original point
				#xnow = xrandom
				#ynow = yrandom
				#print(ncircles)

				
				#We found a boundary, so we break this walk
			#	break
	
	
	Uestimate = 0
	for FoundBpoint in FoundBoundaries:
		Uestimate += FoundBpoint[2]
	Uestimate *= 1/len(FoundBoundaries)
	print(len(FoundBoundaries))
	Upoints.append([xrandom,yrandom,Uestimate])
		

#print("Number of found Upoints")
#print
	



X = np.zeros(len(Upoints))
Y = np.zeros(len(Upoints))
Z = np.zeros(len(Upoints))
for i in range(len(Upoints)):
	X[i] = Upoints[i][0]
	Y[i] = Upoints[i][1]
	Z[i] = Upoints[i][2]
	#print(Upoints[i][0],Upoints[i][1],Upoints[i][2])
	#ax.scatter(Upoints[i][0],Upoints[i][1],Upoints[i][2],color = "blue")

Zi = griddata(X,Y,Z,x,y,interp="linear")
#plt.contourf(x,y,Zi)
#zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

Xm,Ym = np.meshgrid(x,y)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Random Walk On Sphere Method")
ax.plot_surface(Xm, Ym, Zi)#,cmap = 'winter')

plt.show()


##########################################################################
#Do circle plot

#Modify epsilon
epsilon = 0.08

#One single circle, starting from one single point
Circlepoints = []

#Here we pick a point on the grid
xrandom = 0.65#np.random.uniform(0,Lx)
yrandom = 0.65#np.random.uniform(0,Ly)


FoundBoundaries = []#np.array([])
#Either loops or nested function maybe

xnow = xrandom
ynow = yrandom

FoundBoundary = False
while FoundBoundary == False:

	#Måske prøv at plotte dem... like, append disse points til en list, og plot dem
	#For at se om de kommer MEGET FAR AWAY from boundary, eller hvad der lige foregår...
	#Måske er det noget floating point error
	#if xnow > 1 or xnow < 0 or ynow > 1 or ynow < 0:
	#	print("Point is outside boundary!")


	#We check first if current point is near/at boundary
	#FoundBoundary = False
	for Bpoint in Boundary:
		d = np.sqrt((xnow-Bpoint[0])**2+(ynow-Bpoint[1])**2)
		
		if d < eps:
			#We found a boundary, and we can end this for loop and later the walk
			FoundBoundaries.append(Bpoint)
			FoundBoundary = True
			
			Circlepoints.append(Bpoint)
			
			#We're also gonna break the loop, there's no need to go through ALL the rest of the points
			break
	
	if FoundBoundary == False:
		#Find biggest allowable radius, is also the minimum distance to boundary
		dmin = Lx
		for Bpoint in Boundary:
			d = np.sqrt((xnow-Bpoint[0])**2+(ynow-Bpoint[1])**2)
			if d < dmin:
				dmin = d
		
		#We add the point we're at, and the max radius  from this point.
		Circlepoints.append([xnow,ynow, dmin])
		
		#We go to a new point
		#New point is xnew = xold + dx, ynew = yold + dy
		anglerandom = np.random.uniform(0,2*np.pi)
		
		xnow = xnow + dmin*np.cos(anglerandom)
		ynow = ynow + dmin*np.sin(anglerandom)

print(Circlepoints)
fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Walk on spheres")
ax.plot([0,0],[0,Ly],color="black")
ax.plot([Lx,Lx],[0,Ly],color="black")
ax.plot([0,Lx],[0,0],color="black")
ax.plot([0,Lx],[Ly,Ly],color="black")
ax.axis([0-0.1,Lx+0.1,0-0.1,Ly+0.1])

Xc = np.zeros(len(Circlepoints))
Yc = np.zeros(len(Circlepoints))
for i in range(len(Circlepoints)):
	Xc[i] = Circlepoints[i][0]
	Yc[i] = Circlepoints[i][1]
	plt.scatter(Circlepoints[i][0],Circlepoints[i][1],color="black")
	
	if i != len(Circlepoints)-1:
		circle1 = plt.Circle((Circlepoints[i][0], Circlepoints[i][1]), Circlepoints[i][2],color='r',fill=False)
		ax.add_artist(circle1)

ax.plot(Xc,Yc,color="black")

plt.show()
	
	
#################################################
#MAKE SUBPLOT
#fig, ax = plt.subplots()
fig = plt.figure(figsize=plt.figaspect(0.5))
#plt.title("Random Walk On Spheres Method")

fig.suptitle("Random Walk On Spheres Method", fontsize=16)

ax = fig.add_subplot(121,projection='3d')

ax.set_xlabel("x")
ax.set_ylabel("y")
#ax.set_title("Random Walk On Sphere Method")
ax.plot_surface(Xm, Ym, Zi)#,cmap = 'winter')



#plt.subplot(122)
ax = fig.add_subplot(122)
ax.plot([0,0],[0,Ly],color="black",alpha=0.5)
ax.plot([Lx,Lx],[0,Ly],color="black",alpha=0.5)
ax.plot([0,Lx],[0,0],color="black",alpha=0.5)
ax.plot([0,Lx],[Ly,Ly],color="black",alpha=0.5)
for Bpoint in Boundary:
	plt.scatter(Bpoint[0],Bpoint[1],color="black")
	
ax.axis([0-0.1,Lx+0.1,0-0.1,Ly+0.1])
ax.set_xlabel("x")
ax.set_ylabel("y")
for i in range(len(Circlepoints)-1):
	#plt.scatter(Circlepoints[i][0],Circlepoints[i][1],color="red")
	
	#ax.arrow(Circlepoints[i][0], Circlepoints[i][1],
	#		Circlepoints[i+1][0]-Circlepoints[i][0], Circlepoints[i+1][1]-Circlepoints[i][1], 
	#		head_width=0.03, head_length=0.05, fc='k', ec='k')
	
	#ax.quiver(Circlepoints[i][0], Circlepoints[i][1],
	#Circlepoints[i+1][0]-Circlepoints[i][0], Circlepoints[i+1][1]-Circlepoints[i][1])
	
	#units='width')
	#circle1 = plt.Circle((Circlepoints[i][0], Circlepoints[i][1]), Circlepoints[i][2],color='r',fill=False)
	#ax.add_artist(circle1)
	
	ax.annotate('', (Circlepoints[i+1][0],Circlepoints[i+1][1]), (Circlepoints[i][0], Circlepoints[i][1]),
				arrowprops=dict(arrowstyle='->'),size=15)
	
	if i != len(Circlepoints)-1:
		circle1 = plt.Circle((Circlepoints[i][0], Circlepoints[i][1]), Circlepoints[i][2],color='red',fill=False)
		ax.add_artist(circle1)

	
plt.scatter(Circlepoints[-1][0],Circlepoints[-1][1],color="red")
ax.plot(Xc,Yc,color="black")#,marker="v")


fig.savefig('WalkOnSpheresSubplot.png', bbox_inches='tight')
plt.show()