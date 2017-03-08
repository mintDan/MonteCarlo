import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


nx = 75
ny = 75
#x = np.linspace(0,10,nx)
#y = np.linspace(0,10,ny)
T = 0.1
J = 1# J > 0 gives ferromagnetism

indices = [i for i in range(nx)]

#X,Y = np.meshgrid(x,y)

#Z = np.random.random_integers(-1,1, size=(10, 10))

S = np.zeros((ny,nx),dtype=int)+1

def randomS(S):
	#Start with random initial spin
	for i in range(nx):
		for j in range(ny):
			val = np.random.rand()
			if val > 0.5:
				S[j,i] = 1
			else:
				S[j,i] = -1
randomS(S)
print(S)

M = np.sum(S)


def CalcH(S):
	
	H = 0

	#Jeg bør bruge nogle pointers, så jeg ikke double count.. fx "S12*S11" eksisterer, så "S11*S12" skal ikke 
	#tilføjes til H... Så måske først lave pairs, og derefter sum dem up..
	#Eller, divide by 4
	#Pairs = []
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


	H = H/2 #How many times does it count the same? Gotta try it 4x4 grid maybe
	#Energien burde være integer så vidt jeg kan se, men lige nu får jeg 190.5 etc... de der 0.5
	#Kan jo prøve med grid i full +1, så bør den give en easy analytical expression for energy som
	#jeg kan compare...
	#Jeg tror dog højest det bliver double counted, ikke 4x counted...
	return H

	
def Esite(S,j,i):
	
	E1 = 0
	
	if j == 0:
		E1 += -J*S[nx-1,i]*S[j,i]
	elif j == nx-1:
		E1 += -J*S[0,i]*S[j,i]
	else:
		E1 += -J*S[j-1,i]*S[j,i]
		E1 += -J*S[j+1,i]*S[j,i]
		

	if i == 0:
		E1 += -J*S[j,ny-1]*S[j,i]
	elif i == ny-1:
		E1 += -J*S[j,0]*S[j,i]
	else:
		E1 += -J*S[j,i+1]*S[j,i]
		E1 += -J*S[j,i-1]*S[j,i]
	
	return 2*E1

def Esiteflip(S,j,i):
	
	E2 = 0
	
	#Outcommented J, we work with natural units anyway
	#Should see if this energy is correct
	
	if j == 0:
		E2 += -S[nx-1,i]*(-1)
		E2 += -S[1,i]*(-1)
	elif j == nx-1:
		E2 += -S[0,i]*(-1)
		E2 += -S[nx-2,i]*(-1)
	else:
		E2 += -S[j-1,i]*(-1)
		E2 += -S[j+1,i]*(-1)
		

	if i == 0:
		E2 += -S[j,ny-1]*(-1)
		E2 += -S[j,1]*(-1)
	elif i == ny-1:
		E2 += -S[j,0]*(-1)
		E2 += -S[j,ny-2]*(-1)
	else:
		E2 += -S[j,i+1]*(-1)
		E2 += -S[j,i-1]*(-1)
	
	#Da summen bliver J*(term+term+term+term), så får vi 4J, men vi kan gå fra -8 til 8, right?
	#Så jeg skal calculate 16 exponentials? Eller 17, med 0...
	#Så hvis vi siger at PreCalcExp[8] = np.exp(0
	#Eller er det noget med at dE altid er et lige tal? Dette vil gøre tingene mere simple
	#Ja, det er jo altid lige tal, forid vi har faktor 2... så jeg kan fjerne nogle af precalc
	#exponentials, men, whatever...
	return 2*E2*S[j,i]
	

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



#Precalculated exponentials
PreCalcExp = [np.exp(-(i-8.0)/T) for i in range(17)]
#exp(-(0-8)/T)
#exp(-(1-8)/T)
#exp(-(2-8)/T)
#exp(-(3-8)/T)
#exp(-(4-8)/T)
#exp(-(5-8)/T)
#exp(-(6-8)/T)
#exp(-(7-8)/T)
#exp(-(8-8)/T)
#exp(-(9-8)/T)
#exp(-(10-8)/T)
#exp(-(11-8)/T)
#exp(-(12-8)/T)
#exp(-(13-8)/T)
#exp(-(14-8)/T)
#exp(-(15-8)/T)
#exp(-(16-8)/T)
#print(PreCalcExp[0])
#print(PreCalcExp)



fig = plt.figure()
ax = fig.gca()
plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
ax.set_title('Ising model, Metropolis algorithm')
fig.savefig('Ising.png', bbox_inches='tight')
#plt.colorbar()
def animate(i): #i increment with 1 each step

	ax.clear()
	ax.set_title('Ising model, Metropolis algorithm')
	
	#Vi skal update den gamle energy til den nye energy, right? som x=xold
	#H1 = CalcH(S)
	#M = np.sum(S)
	#print("E = {} M = {}".format(H1,M))
	
	
	randi = np.random.randint(0,nx)
	randj = np.random.randint(0,nx)
	#print("row,column")
	#print(randj,randi)
	#print(S)
	#print("Next")
	#E1 = Esite(S,randi,randj)
	E2 = Esiteflip(S,randj,randi)
	#print(E2)
	
	#Stry = S.copy()
	#Stry[randj,randi] = -1*Stry[randj,randi]
	#print(Stry)

	#H2 = CalcH(Stry)

	#dH = H2-H1
	#print("Change in energy")
	#print(dH)
	
	dH = E2#-E1
	
	if dH < 0:
		S[randj,randi] = -1*S[randj,randi]
	else:
		x = np.random.uniform(0,1)
		
		P = PreCalcExp[dH-8]#
		#Hvis dH=0, så tager vi PreCalcExp[8]
		#P = np.exp(-dH/T)
		if x <= P:
			S[randj,randi] = -1*S[randj,randi]

	
	


	plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
	#plt.colorbar()
	return None
anim = animation.FuncAnimation(fig, animate, frames=10, interval=100)

plt.show()
####################################
#jeg kan prøve at se med 3x3 case og se om den calculate det samme som jeg gør på paper, og måske
#4x4


