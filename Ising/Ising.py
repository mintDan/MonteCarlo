import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


nx = 20
ny = 20
x = np.linspace(0,10,nx)
y = np.linspace(0,10,ny)
T = 0.3

indices = [i for i in range(nx)]

X,Y = np.meshgrid(x,y)

#Z = np.random.random_integers(-1,1, size=(10, 10))

S = np.zeros((ny,nx))+1

def randomS(S):
	for i in range(nx):
		for j in range(ny):
			val = np.random.rand()
			if val > 0.1:
				S[j,i] = 1
			else:
				S[j,i] = -1
randomS(S)
print(S)

M = np.sum(S)


def CalcH(S):
	J = 1 # J > 0 gives ferromagnetism
	H = 0

	#Jeg bør bruge nogle pointers, så jeg ikke double count.. fx "S12*S11" eksisterer, så "S11*S12" skal ikke 
	#tilføjes til H... Så måske først lave pairs, og derefter sum dem up..
	#Eller, divide by 4
	Pairs = []
	for i in range(nx):
		for j in range(ny):

			
			#for n in range(-1,2):
			#	for m in range(-1,2):

					#Måske er der en bedre måde at exclude en value fra range?
					#if n != 0 and m != 0:
					#	if j+n in indices and i+m in indices:

					#		Spin1 = "S{}{}".format(j+n,i+m)
					#		Spin2 = "S{}{}".format(j,i)
							
					#		pair1 = [Spin1,Spin2]
					#		pair2 = [Spin2,Spin1]

					#		if pair1 not in Pairs and pair2 not in Pairs:
					#			Pairs.append(pair1)
					#			Pairs.append(pair2)
					#			H += -J*S[j+n,i+m]*S[j,i]




			#Nearest neighbor? Ingen diagonal neighbors?
			#PDE, periodic boundary conditions
			if j == 0:
				H += -J*S[nx-1,i]*S[j,i]
			elif j == nx-1:
				H += -J*S[0,i]*S[j,i]
			else:
				H += -J*S[j-1,i]*S[j,i]
				H += -J*S[j+1,i]*S[j,i]
				
			
			if i == 0:
				H += -J*S[j,ny-1]*S[j,i]
			elif i == ny-1:
				H += -J*S[j,0]*S[j,i]
			else:
				H += -J*S[j,i+1]*S[j,i]
				H += -J*S[j,i-1]*S[j,i]
			
			

	#Noget med at jeg skal calculate probabilities, right? Eller, nej, måske ikke, men man KUNNE godt...
	#Men der vil være 10^10 mulige permutations så hvis jeg tænker rigtigt...
	#Hvis man vil brute force probabilities...

	#print(Pairs)
	#print(H)
	#print(M)
	H = H/4
	return H
H1 = CalcH(S)

#Anyway, nu skal vi flip spin af en random thingy right?
randi = np.random.randint(0,nx)
randj = np.random.randint(0,nx)
print("row,column")
print(randj,randi)
print(S)
print("Next")
Stry = S.copy()
Stry[randj,randi] = -1*Stry[randj,randi]
print(Stry)

H2 = CalcH(Stry)

dH = H2-H1
print("Change in energy")
print(dH)

if dH < 0:
	S[randj,randi] = -1*S[randj,randi]
else:
	x = np.random.uniform(0,1)
	P = np.exp(-dH/T)
	if x <= P:
		S[randj,randi] = -1*S[randj,randi]




fig = plt.figure()
ax = fig.gca()
plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
ax.set_title('Ising model, Metropolis algorithm')
fig.savefig('Ising.png', bbox_inches='tight')
#plt.colorbar()
def animate(i): #i increment with 1 each step

	ax.clear()
	ax.set_title('Ising model, Metropolis algorithm')

	H1 = CalcH(S)
	print(H1)
	randi = np.random.randint(0,nx)
	randj = np.random.randint(0,nx)
	#print("row,column")
	#print(randj,randi)
	#print(S)
	#print("Next")
	Stry = S.copy()
	Stry[randj,randi] = -1*Stry[randj,randi]
	#print(Stry)

	H2 = CalcH(Stry)

	dH = H2-H1
	#print("Change in energy")
	#print(dH)

	if dH < 0:
		S[randj,randi] = -1*S[randj,randi]
	else:
		x = np.random.uniform(0,1)
		P = np.exp(-dH/T)
		if x <= P:
			S[randj,randi] = -1*S[randj,randi]

	#Vi skal update den gamle energy til den nye energy, right? som x=xold
	


	plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
	#plt.colorbar()
	return None
anim = animation.FuncAnimation(fig, animate, frames=10, interval=100)

plt.show()
####################################
#jeg kan prøve at se med 3x3 case og se om den calculate det samme som jeg gør på paper, og måske
#4x4


