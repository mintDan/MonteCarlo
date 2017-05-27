import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import multiprocessing
import time


def TimeCode(f):
	"""
	Decorator to time code
	"""
	def timed(*args, **kw):
		ts = time.time()
		result = f(*args, **kw)
		te = time.time()
		print("Time difference based on decorator")
		print(f.__name__)
		print(te-ts)
		
		return result #Den return result fra RunMonteCarlo()
	
	#Denne her return en function, ikke et tal etc, men en function, important difference
	return timed #Den return result fra timed, som return result fra RunMonteCarlo()

def CalcHSlow(S,nx,ny,J):
	"""
	Calculates Hamiltonian, but this is inefficient. Rather the function Esiteflip(S,i,j) is used
	"""
	
	H = 0
	
	#H2 = 0

	#Jeg bør bruge nogle pointers, så jeg ikke double count.. fx "S12*S11" eksisterer, så "S11*S12" skal ikke 
	#tilføjes til H... Så måske først lave pairs, og derefter sum dem up..
	#Eller, divide by 4
	#Pairs = []
	#Jeg udkommenterer de nederste, og lægger dem til de øverste. Færre operations tror jeg.
	
	for i in range(nx):
		for j in range(ny):

			#Nearest neighbor? Ingen diagonal neighbors?
			#PDE, periodic boundary conditions
			if j == 0:
				H += S[nx-1,i]*S[j,i]+S[1,i]*S[j,i]
				#H += S[1,i]*S[j,i]
				
				#H2 += np.abs(-J*S[nx-1,i]*S[j,i])
				#H2 += np.abs(-J*S[1,i]*S[j,i])
			elif j == nx-1:
				H += S[0,i]*S[j,i]+S[nx-2,i]*S[j,i]
				#H += S[nx-2,i]*S[j,i]
				
				#H2 += np.abs(-J*S[0,i]*S[j,i])
				#H2 += np.abs(-J*S[nx-2,i]*S[j,i])
			else:
				H += S[j-1,i]*S[j,i]+S[j+1,i]*S[j,i]
				#H += S[j+1,i]*S[j,i]
				
				#H2 += np.abs(-J*S[j-1,i]*S[j,i])
				#H2 += np.abs(-J*S[j+1,i]*S[j,i])
				
			
			if i == 0:
				H += S[j,ny-1]*S[j,i]+S[j,1]*S[j,i]
				#H += S[j,1]*S[j,i]
				
				#H2 += np.abs(-J*S[j,ny-1]*S[j,i])
				#H2 += np.abs(-J*S[j,1]*S[j,i])
				
			elif i == ny-1:
				H += S[j,0]*S[j,i]+S[j,ny-2]*S[j,i]
				#H += S[j,ny-2]*S[j,i]
				
				#H2 += np.abs(-J*S[j,0]*S[j,i])
				#H2 += np.abs(-J*S[j,ny-2]*S[j,i])
			else:
				H += S[j,i+1]*S[j,i]+S[j,i-1]*S[j,i]
				#H += S[j,i-1]*S[j,i]
				
				#H2 += np.abs(-J*S[j,i+1]*S[j,i])
				#H2 += np.abs(-J*S[j,i-1]*S[j,i])
	
	#Noget med at jeg skal calculate probabilities, right? Eller, nej, måske ikke, men man KUNNE godt...
	#Men der vil være 10^10 mulige permutations så hvis jeg tænker rigtigt...
	#Hvis man vil brute force probabilities...

	
	H = -J*H/2

	#H2 = H2/2 #How many times does it count the same? Gotta try it 4x4 grid maybe
	#Energien burde være integer så vidt jeg kan se, men lige nu får jeg 190.5 etc... de der 0.5
	#Kan jo prøve med grid i full +1, så bør den give en easy analytical expression for energy som
	#jeg kan compare...
	#Jeg tror dog højest det bliver double counted, ikke 4x counted...
	return H
	
def CalcH(S,nx,ny,J):
	"""
	Calculates Hamiltonian, but this is inefficient. Rather the function Esiteflip(S,i,j) is used
	"""
	
	H = 0
	
	#H2 = 0

	#Jeg bør bruge nogle pointers, så jeg ikke double count.. fx "S12*S11" eksisterer, så "S11*S12" skal ikke 
	#tilføjes til H... Så måske først lave pairs, og derefter sum dem up..
	#Eller, divide by 4
	#Pairs = []
	
	#Jeg udkommenterer de nederste og lægger dem direkte til de øverste, færre operations
	for i in range(nx):
		for j in range(i+1,ny):

			#Nearest neighbor? Ingen diagonal neighbors?
			#PDE, periodic boundary conditions
			if j == 0:
				H += S[nx-1,i]*S[j,i]+S[1,i]*S[j,i]
				#H += S[1,i]*S[j,i]
				
				#H2 += np.abs(-J*S[nx-1,i]*S[j,i])
				#H2 += np.abs(-J*S[1,i]*S[j,i])
			elif j == nx-1:
				H += S[0,i]*S[j,i]+S[nx-2,i]*S[j,i]
				#H += S[nx-2,i]*S[j,i]
				
				#H2 += np.abs(-J*S[0,i]*S[j,i])
				#H2 += np.abs(-J*S[nx-2,i]*S[j,i])
			else:
				H += S[j-1,i]*S[j,i]+S[j+1,i]*S[j,i]
				#H += S[j+1,i]*S[j,i]
				
				#H2 += np.abs(-J*S[j-1,i]*S[j,i])
				#H2 += np.abs(-J*S[j+1,i]*S[j,i])
				
			
			if i == 0:
				H += S[j,ny-1]*S[j,i]+S[j,1]*S[j,i]
				#H += S[j,1]*S[j,i]
				
				#H2 += np.abs(-J*S[j,ny-1]*S[j,i])
				#H2 += np.abs(-J*S[j,1]*S[j,i])
				
			elif i == ny-1:
				H += S[j,0]*S[j,i]+S[j,ny-2]*S[j,i]
				#H += S[j,ny-2]*S[j,i]
				
				#H2 += np.abs(-J*S[j,0]*S[j,i])
				#H2 += np.abs(-J*S[j,ny-2]*S[j,i])
			else:
				H += S[j,i+1]*S[j,i]+S[j,i-1]*S[j,i]
				#H += S[j,i-1]*S[j,i]
				
				#H2 += np.abs(-J*S[j,i+1]*S[j,i])
				#H2 += np.abs(-J*S[j,i-1]*S[j,i])
	
	#Noget med at jeg skal calculate probabilities, right? Eller, nej, måske ikke, men man KUNNE godt...
	#Men der vil være 10^10 mulige permutations så hvis jeg tænker rigtigt...
	#Hvis man vil brute force probabilities...

	
	H = -J*H

	#H2 = H2/2 #How many times does it count the same? Gotta try it 4x4 grid maybe
	#Energien burde være integer så vidt jeg kan se, men lige nu får jeg 190.5 etc... de der 0.5
	#Kan jo prøve med grid i full +1, så bør den give en easy analytical expression for energy som
	#jeg kan compare...
	#Jeg tror dog højest det bliver double counted, ikke 4x counted...
	return H

def Esiteflip(S,j,i,nx,ny):
	"""
	More efficient way of calculating dH
	If this can't be used with Multiprocessing then just deifne this function inside the other function
	"""
	
	E2 = 0
	
	#Outcommented J, we work with natural units anyway
	#Should see if this energy is correct
	
	
	#Jeg udkommenterer de nederste og lægger dem til de øverste, laver nok lidt mindre operations på den måde.
	if j == 0:
		E2 += -S[nx-1,i]-S[1,i]#*(-1)
		#E2 += -S[1,i]#*(-1)
	elif j == nx-1:
		E2 += -S[0,i]-S[nx-2,i]#*(-1)
		#E2 += -S[nx-2,i]#*(-1)
	else:
		E2 += -S[j-1,i]-S[j+1,i]#*(-1)
		#E2 += -S[j+1,i]#*(-1)
		

	if i == 0:
		E2 += -S[j,ny-1]-S[j,1]#*(-1)
		#E2 += -S[j,1]#*(-1)
	elif i == ny-1:
		E2 += -S[j,0]-S[j,ny-2]#*(-1)
		#E2 += -S[j,ny-2]#*(-1)
	else:
		E2 += -S[j,i+1]-S[j,i-1]#*(-1)
		#E2 += -S[j,i-1]#*(-1)
	
	#Da summen bliver J*(term+term+term+term), så får vi 4J, men vi kan gå fra -8 til 8, right?
	#Så jeg skal calculate 16 exponentials? Eller 17, med 0...
	#Så hvis vi siger at PreCalcExp[8] = np.exp(0
	#Eller er det noget med at dE altid er et lige tal? Dette vil gøre tingene mere simple
	#Ja, det er jo altid lige tal, forid vi har faktor 2... så jeg kan fjerne nogle af precalc
	#exponentials, men, whatever...
	return 2*E2*S[j,i]*(-1)
	

# H1 = CalcH(S)

# #Anyway, nu skal vi flip spin af en random thingy right?
# randi = np.random.randint(0,nx)
# randj = np.random.randint(0,nx)
# print("row,column")
# print(randj,randi)
# print(S)
# print("Next")
# Stry = S.copy()
# Stry[randj,randi] = -1*Stry[randj,randi]
# print(Stry)

# H2 = CalcH(Stry)

# dH = H2-H1
# print("Change in energy")
# print(dH)

# if dH < 0:
	# S[randj,randi] *= -1
# else:
	# x = np.random.uniform(0,1)
	# P = np.exp(-dH/T)
	# if x <= P:
		# S[randj,randi] *= -1



#Precalculated exponentials
#PreCalcExp = [np.exp(-(i-8.0)/T) for i in range(17)]
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

#P = PreCalcExp[int(dE)+8]


#PreCalcExp = np.exp([-(2*i-8.0)/T for i in range(9)])
#PreCalcExp[0] = exp(-(0-8)/T)
#PreCalcExp[1] = exp(-(2-8)/T)
#PreCalcExp[2] = exp(-(4-8)/T)
#PreCalcExp[3] = exp(-(6-8)/T)
#PreCalcExp[4] = exp(-(8-8)/T)
#PreCalcExp[5] = exp(-(10-8)/T)
#PreCalcExp[6] = exp(-(12-8)/T)
#PreCalcExp[7] = exp(-(14-8)/T)
#PreCalcExp[8] = exp(-(16-8)/T)



def MCrunParallel(T,N,J,nx,ny,PreCalcExp,out_list):
	"""
	This can be called with multiprocessing function.
	But, doing autocorrelation seems to be troublesome...
	We got 4 arrays that will be copied back and forth, instead of just 4 values...
	So maybe for autocorrelation, make alternative function/method/do it without multiprocessing
	"""

	
	#==================================================================
	#Think i should seed each run? WIthout any input it will use current time to seed.
	np.random.seed()
	
	#==================================================================
	#Here I make the initial spin state. If T < Tcritical, i make fully aligned state of ones, else i make a random state of {-1,1}
	
	#randomS(S)
	if T < 2.2:
		S = np.ones((ny,nx))


	else:
		#If val of index = 0, then we go to -1, if val of index = 2, then we get 1.
		S = 2*np.random.randint(0,2,(ny,nx))-1
		
		
	#print(S)
	
	#Lige her bør jeg faktisk lave en endnu en loop
	
	#Der bør være en Meqold = 0 og 1, 0 for T > Tcrit, det er mere efficient, så vi tager Meqold som givet ud fra initial state her
	#Det er mere flexible
	#Meqold = 1
	Meqold = np.abs(np.sum(S))
	#if Meqold == 0:
		#Jeg må ikke have division by 0 error
	#	Meqold = 0.1
	
	
	MEquilibriumCount = 0
	#100*N = 100*625 = 62500
	#E2calcavg = 0
	

	
	
	#I need to make an average of M*M for <M**2>
	#And later, I need to make an average of THESE aswell, with 1/ntest faktor
	#So it will be double average for <M**2>, and probably should do the same for <E**2>?
	Avgcount = 0
	
	MavgMC = 0
	M2avgMC = 0
	EavgMC = 0
	E2avgMC = 0
	
	
	#If we get 10 equilibrium points in a row, these will be added to MavgMC etc...
	#Otherwise, if we get say 3 equilibrium points, then 4th test is a failure, then we reset these
	#that way we don't add values from different "equilibrium extremums"
	#MavgMCtemporary = 0
	#M2avgMCtemporary = 0
	#EavgMCtemporary = 0
	#E2avgMCtemporary = 0
	
	
	#===============================================
	#Moving average, using n = 10 points
	nmoving = 10
	Mmoving = np.zeros(nmoving)
	Sweeps = 0
	
	Mmovingavg = 0
	
	
	#================================================
	#This is the actual MonteCarlo loop, changing the configuration based on probabilities
	#This should perhaps be a while loop instead... while Avgcount < 10
	n = 0
	while Avgcount < 10:
	#for n in range(150*N):
		randi = np.random.randint(0,nx)
		randj = np.random.randint(0,ny)


		dE = Esiteflip(S,randj,randi,nx,ny)
	
		if dE < 0:
			S[randj,randi] *= -1
		else:
			x = np.random.uniform(0,1)
			
			#P = np.exp(-dE/T) #PreCalcExp[dH-8]# #Cant use precalculated that well, since we are increasing T
			#P = PreCalcExp[int(dE)+8]
			P = PreCalcExp[int((dE+8)/2)]
			#Hvis dH=0, så tager vi PreCalcExp[8]
			#P = np.exp(-dH/T)
			if x <= P:
				S[randj,randi] *= -1
				
		#===========================================================
		#Autocorrelation
				
		# #Method 2
		# if nt == Nt-1:
			# if navg == ntest-1:
				# if n <= nautotimes-1:
					# Mauto = np.abs(np.sum(S))
					# Marray[n] = Mauto
					# M2array[n] = Mauto*Mauto
					
					# Eauto = CalcH(S)
					# Earray[n] = Eauto
					# E2array[n] = Eauto*Eauto
					
		#=========================================================
		#Here I count succesive equilibrium states,
		#It needs to "remember" at least 2 magnetizations, the current and the last
		#Can be made more advanced, but that's the simplest
		#Samples at N,2N,3N,4N... monte carlo steps
		if n%N == 0:
			#We only measure M during full Monte Carlo Sweeps, every N steps iirc
			
			Sweeps += 1
			#Meq for Mequilibriation
			Meq = np.abs(np.sum(S))
			
			
			#Let's make moving average.
			#The moving average is just moving average.
			#It doesn't reset etc like the equilibriumcounters, it's always there.
			

			#Sweeps = 1 -> Mmoving[0] = Meq
			#Sweeps = 2 -> Mmoving[1] = Meq
			#Sweeps = 3 -> Mmoving[2] = Meq
			#Sweeps = 4 -> Mmoving[3] = Meq
			#Sweeps = 5 -> Mmoving[4] = Meq
			#Sweeps = 6 -> Mmoving[0] = Meq
			#Sweeps = 7 -> Mmoving[1] = Meq
			#Sweeps = 8 -> Mmoving[2] = Meq
			#Sweeps = 9 -> Mmoving[3] = Meq
			#Sweeps = 10 -> Mmoving[4] = Meq
			#Sweeps = 11 -> Mmoving[0] = Meq
			
			#nmoving = 5
			#6/5 = 1 + 1/5
			#7/5 = 1 + 2/5
			
			movingindex = (Sweeps%nmoving)-1
			
			Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[movingindex]/nmoving
			
			#Skal lige finde ud af det med moving average lidt mere...
			
			#Sweeps = 1 -> Mmoving[0] = Mmovingavg + Meq/nmoving - Mmoving[4]
			#Sweeps = 2 -> Mmoving[1] = Mmovingavg + Meq/nmoving - Mmoving[3]
			#Sweeps = 3 -> Mmoving[2] = Mmovingavg + Meq/nmoving - Mmoving[2]
			#Sweeps = 4 -> Mmoving[3] = Mmovingavg + Meq/nmoving - Mmoving[1]
			#Sweeps = 5 -> Mmoving[4] = Mmovingavg + Meq/nmoving - Mmoving[0]
			#Sweeps = 6 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[0]
			#Sweeps = 7 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[1]
			#Sweeps = 8 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[2]
			#Sweeps = 9 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[3]
			#Sweeps = 10 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[4]
			#Sweeps = 11 -> Mmovingavg = Mmovingavg + Meq/nmoving - Mmoving[0]
			
			
						
			Mmoving[movingindex] = Meq
			
			if Meqold != 0:
				#Skal jeg måske lige undersøge hvad der sker, if it so happens, at faktisk Meqold = 0? Men det er ret usandsynligt
				
				#Tror denne her skal ændres tbh? Hvis jeg vil sammenligne to ting, så er det (a-b)/b tror jeg?
				#if np.abs(Meq/(Meqbefore)) < 5:
				#Actually, even with < 0.99, this can take significant timesteps.
				#Because, sometimes old vs current can be, say, 500% appart.
				#Especially at high T, with a lot of randomization, you'll have situations where
				#Meqold = 0 (completely neutral magnetization)
				#But, Meq = 5, it just happened such that 5 spin sites were skewed in one direction.
				#Then we are 500% above Meqold. So actually, with higher T, this condition should be more relaxed!
				#so, for now, add in 2*0.99 factor
				
				#if np.abs((Meq-Meqold)/Meqold) < 3*0.99:
				if np.abs((Meq-Mmovingavg)/Mmovingavg) < 3*0.99:
				#The last state is within e.g 5% of the current, so we have an Equilibrium count
					MEquilibriumCount += 1
					

				
				else:
					#Reset equilibrium counter
					#Hvis vi har fx count == 2, men så kommer der en huge spike, så skal count gå tilbage til 0...
					#We should in fact reset both MEquilibriumCount AND Avgcount... Perhaps there is only need for one of these actually
					MEquilibriumCount = 0
					
					Avgcount = 0
					
					
					#MavgMCtemporary = 0
					#M2avgMCtemporary = 0
					#EavgMCtemporary = 0
					#E2avgMCtemporary = 0
					MavgMC = 0
					M2avgMC = 0
					EavgMC = 0
					E2avgMC = 0
					
			
			if MEquilibriumCount >= 10:
				#We have reached equilibrium let's say, since the last 10 equilibriums are within 5%
				#On the hand, having many equilibrium counts is good, on the other, sometimes by chance, a state will hit a spike,
				#even after hitting equilibrium. so if we demand 100 equilibrium counts in a row, there is a good chance that at SOME POINT
				#during those 100 counts, a subsequent state would have slightly too different, and we would have to start all over with 100 counts.
				#So 10 seems enough.
				
				#These should actually be reset if MEquilibriumCount is reset. So then we should actually remove Meq etc.
				#Så, there should perhaps be another separate counter, and only AFTER we have gone through 10 Equilibrium counts,
				#do we add them ALL UP at once...
				
				#Også... det bør faktisk måske være en moving average....
				#Jeg bør lave en moving average, og sammenlinge om Meq er within % af DENNE, er lidt mere cool synes jeg
				
				#For equilibrium counts above 10, we calculate average M2, M etc
				MavgMC += Meq
				M2avgMC += Meq*Meq
				
				Ecalc = CalcH(S,nx,ny,J)
				EavgMC += Ecalc
				E2avgMC += Ecalc*Ecalc
				
				
				#Avgcount used for E and M
				Avgcount += 1
				
			
			#if MEquilibriumCount == 20:
				#This may interfere with our while loop? We have two conditions for exiting the loop?
				#So we have done 10x measurements for M2, E2 etc, and now we're ready to exit this state
			#	break
				
			#else:
			Meqold = Meq

		n += 1
	#if navg == ntest-1:
	#	print(n)
	#print(Avgcount)
	#The first, inner average of M**2
	#This is the average of samples from the same MC run. We get 10 samples from each MC run, after the configuration has reached equilibrium
	MavgMC *= (1/Avgcount)
	M2avgMC *= (1/Avgcount)
	E2avgMC *= (1/Avgcount)
	EavgMC *= (1/Avgcount)
	#MavgMC = MavgMCtemporary/Avgcount
	#M2avgMC = M2avgMCtemporary/Avgcount
	#E2avgMC = E2avgMCtemporary/Avgcount
	#EavgMC = EavgMCtemporary/Avgcount
	
	
	TupleMCvals = (MavgMC,M2avgMC,E2avgMC,EavgMC)
	out_list.put(TupleMCvals)
	#return MavgMC,M2avgMC,E2avgMC,EavgMC
	#return TupleMCvals
			
			
	


#======================================================
#Make (T,M) graph etc. For each temperature T, we add the final averaged out M,E,XT,Cv
#Ms = []
#Es = []
#XTs = []
#Cvs = []
#Ts = []

@TimeCode
def MainScript(T,Ms,Es,XTs,Cvs,Ts,CEts,CMts):
	"""
	Maybe calling multiprocessing here would be even more efficient. If possible. For lower overhead. 
	So calling multiprocessing one step higher.
	Split the temperature ranges into 4, perhaps.
	So each core get a section, say first core get temperatures T1 = [0.7,0.9,1.1,1.3] etc
	"""
	
	#The main outer loop changes the temperature, so Nt is the number of different temperatures we examine
	for nt in range(Nt):

		
		
		#==================================================================
		#These will be averaged out over ntest times...
		#With independent starting points
		Eavg = 0
		Mavg = 0
		M2avg = 0
		E2avg = 0
		
		
		#==================================================================
		#Precalculated exponenstial, for every T, for more efficiency
		#Right now i make 17, but there's only 5 i think! Because of 2* faktor, so only even numbers
		#T changes each timestep, so we also calculate these each timestep
		#PreCalcExp = [np.exp(-(i-8.0)/T) for i in range(17)]
		#PreCalcExp = np.exp([-(i-8.0)/T for i in range(17)])
		
		
		PreCalcExp = [np.exp(-(2*i-8.0)/T) for i in range(9)]
		
		
		#==================================================================
		#For autocorrelation
		#mt0 = 0
		#m2t0 = 0
		#mt = 0
		#mt0mt = 0
		
		
		Marray = np.zeros(nautotimes)
		M2array = np.zeros(nautotimes)
		
		Earray = np.zeros(nautotimes)
		E2array = np.zeros(nautotimes)
		
		#==================================================================
		#Fordi noget med at bad initial conditions kan give local minim..
		#ntest is the number of times i run the simulation at the SAME temperature, to average results
		#this for loop can perhaps be multiprocessed/parallized?
		#Because, each test is independent
		#So, could call 4 tests at once
		
		#Jeg kunne lave en function MCrun()
		#Som så bliver kaldt ntest gange
		
		#Autocorrelation kan måske være lidt svær at have i multiprocess...
		#Fordi den ændrer jo global array
		
		#==========================================
		#Multiprocessing process
		#Denne her kommer til at holde 4 
		# out_list = multiprocessing.Queue()
		
		# jobs = []
		
		# for i in range(4):
		# #out_list = []
			# numb1 = i*int(Nparticles/4)
			# numb2 = numb1 + int(Nparticles/4)
			# process = multiprocessing.Process(target=Htheorem, 
			                              # args=(SpeedStep,halfmass,numb1,numb2,nE,EnergyLinspace,dE,out_list))
		# #print(process)
		# #print(out_list)
			# jobs.append(process)
		
		# for j in jobs:
			# j.start()
		# #print(out_list)

		# # Ensure all of the processes have finished
		# #result = []
		# for j in jobs:
			# #print(out_list)
			# j.join()
			# #print(out_list)
			# #result.append(j.exitcode)
		# #print(result)
		# results = [out_list.get() for j in jobs]
		# #print(results)
		# Harray[nt] = np.sum(results)
		
		

		
		out_list = multiprocessing.Queue()
		
		jobs = []
		
		for i in range(ntest):
			#out_list = []
			
			process = multiprocessing.Process(target=MCrunParallel, 
											args=(T,N,J,nx,ny,PreCalcExp,out_list))
			jobs.append(process)
											
		for j in jobs:
			j.start()


		# # Ensure all of the processes have finished
		result = []
		for j in jobs:
			j.join()

			result.append(j.exitcode)

		results = [out_list.get() for j in jobs]

		#results = [[(Eavg,MavgMC...)],[(Eavg,MavgMC...)],[(Eavg,MavgMC...)],[(Eavg,MavgMC...)]]
		#results = [(Eavg,MavgMC...),(Eavg,MavgMC...),(Eavg,MavgMC...)]
		#Harray[nt] = np.sum(results)
		
				#Calculate Eavg from equilibration samples
		for i in range(ntest):
			Mavg += results[i][0]
			M2avg += results[i][1]
			E2avg += results[i][2]
			Eavg += results[i][3]

		
		
		
		
		# for navg in range(ntest):
			
			# #==================================================================
			# #Think i should seed each run? WIthout any input it will use current time to seed.
			# np.random.seed()
			
			# #==================================================================
			# #Here I make the initial spin state. If T < Tcritical, i make fully aligned state of ones, else i make a random state of {-1,1}
			
			# #randomS(S)
			# if T < 2.2:
				# S = np.ones((ny,nx))

		
			# else:
				# #If val of index = 0, then we go to -1, if val of index = 2, then we get 1.
				# S = 2*np.random.randint(0,2,(ny,nx))-1
				
				
			# #print(S)
			
			# #Lige her bør jeg faktisk lave en endnu en loop
			
			# #Der bør være en Meqold = 0 og 1, 0 for T > Tcrit, det er mere efficient, så vi tager Meqold som givet ud fra initial state her
			# #Det er mere flexible
			# #Meqold = 1
			# Meqold = np.abs(np.sum(S))
			# if Meqold == 0:
				# #Jeg må ikke have division by 0 error
				# Meqold = 0.1
			
			
			# MEquilibriumCount = 0
			# #100*N = 100*625 = 62500
			# #E2calcavg = 0
			

			
			
			# #I need to make an average of M*M for <M**2>
			# #And later, I need to make an average of THESE aswell, with 1/ntest faktor
			# #So it will be double average for <M**2>, and probably should do the same for <E**2>?
			# Avgcount = 0
			
			# MavgMC = 0
			# M2avgMC = 0
			# EavgMC = 0
			# E2avgMC = 0
			
			# #================================================
			# #This is the actual MonteCarlo loop, changing the configuration based on probabilities
			# #This should perhaps be a while loop instead... while Avgcount < 10
			# n = 0
			# while Avgcount < 10:
			# #for n in range(150*N):
				# randi = np.random.randint(0,nx)
				# randj = np.random.randint(0,nx)


				# dE = Esiteflip(S,randj,randi)
			
				# if dE < 0:
					# S[randj,randi] *= -1
				# else:
					# x = np.random.uniform(0,1)
					
					# #P = np.exp(-dE/T) #PreCalcExp[dH-8]# #Cant use precalculated that well, since we are increasing T
					# #P = PreCalcExp[int(dE)+8]
					# P = PreCalcExp[int((dE+8)/2)]
					# #Hvis dH=0, så tager vi PreCalcExp[8]
					# #P = np.exp(-dH/T)
					# if x <= P:
						# S[randj,randi] *= -1
						
				# #=========================================================
				# #Values for autocorrelation
				# #if nt == Nt-1:
					# #Do autocorrelation
					# #I need values at t0, let's do it for magnetization made
					# #let's say it's at timestep 200 i call it t0
					# #So, the time is an n index and obviously has to be somewhere in the MONTECARLO LOOP.
					# #So, i actually do make a counter for n in the monte carlo loop, so
				# #if n == 200:
				# #	mt0current = np.abs(np.sum(S))
					
				# #	mt0 += mt0current
				# #	m2t0 += mt0current*mt0current
				
				# #if n == 300:
				# #	mtcurrent = np.abs(np.sum(S))
				# #	mt += mtcurrent
					
				# #	mt0mt += mt0current*mtcurrent
					
				
				# #Method 2
				# if nt == Nt-1:
					# if navg == ntest-1:
						# if n <= nautotimes-1:
							# Mauto = np.abs(np.sum(S))
							# Marray[n] = Mauto
							# M2array[n] = Mauto*Mauto
							
							# Eauto = CalcH(S)
							# Earray[n] = Eauto
							# E2array[n] = Eauto*Eauto
						
					
				
				
				# #=========================================================
				# #Here I count succesive equilibrium states,
				# #It needs to "remember" at least 2 magnetizations, the current and the last
				# #Can be made more advanced, but that's the simplest
				# #Samples at N,2N,3N,4N... monte carlo steps
				# if n%N == 0:
				
					# #We only measure M during full Monte Carlo Sweeps, every N steps iirc
					
					# Meq = np.abs(np.sum(S))
					# if Meqold != 0:
						
						# #Tror denne her skal ændres tbh? Hvis jeg vil sammenligne to ting, så er det (a-b)/b tror jeg?
						# #if np.abs(Meq/(Meqbefore)) < 5:
						# if np.abs((Meq-Meqold)/Meqold) < 0.99:
						# #The last state is within 5% of the current, so we have an Equilibrium count
							# MEquilibriumCount += 1
							
						# #Reset equilibrium counter
						# #Hvis vi har fx count == 2, men så kommer der en huge spike, så skal count gå tilbage til 0...
						# else:
							# MEquilibriumCount = 0
					
					# if MEquilibriumCount >= 10:
						# #We have reached equilibrium let's say, since the last 10 equilibriums are within 5%
						# #On the hand, having many equilibrium counts is good, on the other, sometimes by chance, a state will hit a spike,
						# #even after hitting equilibrium. so if we demand 100 equilibrium counts in a row, there is a good chance that at SOME POINT
						# #during those 100 counts, a subsequent state would have slightly too different, and we would have to start all over with 100 counts.
						# #So 10 seems enough.
						
						# #For equilibrium counts above 10, we calculate average M2
						# MavgMC += Meq
						# M2avgMC += Meq*Meq
						
						# Ecalc = CalcH(S)
						# EavgMC += Ecalc
						# E2avgMC += Ecalc*Ecalc
						
						
						# #Avgcount used for E and M
						# Avgcount += 1
						
					
					# if MEquilibriumCount == 20:
					
						# #So we have done 10x measurements for M2, E2 etc, and now we're ready to exit this state
						# break
						
					# else:
						# Meqold = Meq

				# n += 1
			# if navg == ntest-1:
				# print(n)	
				
			# #The first, inner average of M**2
			# #This is the average of samples from the same MC run. We get 10 samples from each MC run, after the configuration has reached equilibrium
			# MavgMC *= (1/Avgcount)
			# M2avgMC *= (1/Avgcount)
			# E2avgMC *= (1/Avgcount)
			# EavgMC *= (1/Avgcount)
					
					
					
			# #E2calcavg *= (1.0/20.0)

			# #MavgMC, M2avgMC, E2avgMC, EavgMC = MCrunParallel(T,nx,ny,PreCalcExp)
			# #==================================================================
			# #What happens here?
			# #Calculate Energies, jeg venter med at divide med N til sidst, pga E2avg
			# #Etotal = CalcH(S)
			
			# #Eavg += Etotal
			
			# #Calculate Eavg from equilibration samples
			# Eavg += EavgMC
			
			
			# #E2avg += Etotal*Etotal
			
			# #Calculate E2avg with the already averaged E2avgMC, to get <E**2>
			# E2avg += E2avgMC
			
			
			# #Calculate Magnetization
			# #I want it to be absolute value, and an average, pr site, i think.
			# #I want to to be absolute especially if I'm gonna run the test multiple times and average it.
			
			# #Mtotal = np.abs(np.sum(S))
			
			# #Mavg += Mtotal
			
			# Mavg += MavgMC
			
			# #Used for isothermal susceptibility, skal måske her divide med N**2?
			# #Skal måske også lave den til M*M ligesom Etotal*Etotal...
			# #I (<M2>-<M>**2) så vil <M>**2 have faktor 1/N**2, mens <M2> kun har faktor 1/N
			# #Så det skal nok være faktor 1/N**2 også på M2
			
			# #M2avg += Mtotal*Mtotal
			
			# #Calculate M2avg with the already averaged M2avgMC, to REALLY get a <M**2>
			
			# M2avg += M2avgMC
			
			# #M2avg = M2avg + np.abs(np.sum(S**2))
		
		
		
		#================================================
		#The energy E and magnetizaion M, etc, for this GIVEN temperature T is calculated.
		#We have done e.g ntest = 10 MC runs at the same temperature,
		#and at each MC run we did e.g Avgcount = 10 succesive samples after reaching equilibrium
		E = Eavg/ntest
		E2 = E2avg/ntest
		
		M = Mavg/ntest
		M2 = M2avg/ntest
		
		
		#Calculate Isothermal susceptibility
		#Det er egentlig XT pr site... men det er pr site squared, virker ikke helt rigtigt
		XT = (1/T)*(M2-M**2)/N
		
		#Calculate specific heat
		#Kan også prøve at sammenligne med finite differences
		#Dette er egentlig Cv pr mass/site
		#Hvilke values af Cv får de andre?
		
		#Lad os køre specific heat uden N... det er jo <E**2>-<E>**2, men IKKE PR SPIN.
		Cv = (1/T**2)*(E2-E**2)/N
		#Det kan godt være, at i stedet for np.abs(), skal det være ( )**2
		
		
		#==============================
		#Energy and magnetization pr spinsite
		#After having used E,M for XT and Cv, now we can divide by N
		E = E/N
		M = M/N
			
		
		
		#Save figure
		if SaveFig == 1:
			fig = plt.figure()
			ax = fig.gca()
			plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
			ax.set_title('Ising model, Metropolis algorithm')
			fig.savefig('Ising{0:0.3}.png'.format(T), bbox_inches='tight')
			plt.close()
		#=================================================================
		#print some stuff to see progress
		print("T = {0:.3f} M = {1:.3f} E = {2:.3f}".format(T,M,E))
		
		#================================================================
		#Append to lists for plotting later
		Ms.append(M)
		Es.append(E)
		XTs.append(XT)
		Ts.append(T)
		Cvs.append(Cv)
		T += dT
		
		
		#mt0 = 0
		#m2t0 = 0
		#mt = 0
		#mt0mt = 0
		
		#==========================================
		#Let's try do see some Autocorrelation, first for last timestep
		#for ntest in range(Ntest):
		#	for nmc in range(BigNmc):
		#		if nt == Nt-1:
					#Do autocorrelation
					#I need values at t0, let's do it for magnetization made
					#let's say it's at timestep 200 i call it t0
					#So, the time is an n index and obviously has to be somewhere in the MONTECARLO LOOP.
					#So, i actually do make a counter for n in the monte carlo loop, so
		#			if n == 200:
		#				mt0current = np.abs(np.sum(S))
		#				
		#				mt0 += mt0current
		#				m2t0 += mt0current*mt0current
					
		#			if n == 300:
		#				mtcurrent = np.abs(np.sum(S))
		#				mt += mtcurrent
						
		#				mt0mt += mt0current*mtcurrent
					
		#mt0avg = mt0/ntest
		#m2t0avg = m2t0/ntest
		#mtavg = mt/ntest
		
		#mt0mtavg = mt0mt/ntest
		
		#Cmt = (mt0mtavg-mt0avg*mtavg)/(m2t0avg-mt0avg*mt0avg)
		#print(Cmt)
		
		
		#Cmt2 = (mt0mtavg-mt0avg*mt0avg)/(m2t0avg-mt0avg*mt0avg)
		#print(Cmt, Cmt2)
		
		
		
		#Method 2
		
		#Marrayavg = np.average(Marray)
		#M2arrayavg = np.average(M2array)
		#Weirdsum = 0
		#for i in range(700):
		#	Weirdsum += Marray[i]*Marray[i+300]
		#Weirdsum *= (1.0/(700))
		
		#C300Method2 = (Weirdsum - Marrayavg*Marrayavg)/(M2arrayavg-Marrayavg*Marrayavg)
		#print(C300)
		
		
		#method 3
		#Weirdsum2 = 0
		#Weirdsum3 = 0
		#for i in range(700):
		#	Weirdsum2 += Marray[i]
		#	Weirdsum3 += Marray[i+300]
		#Weirdsum2 *= (1.0/(700))
		#Weirdsum3 *= (1.0/(700))

		#C300Method3 = Weirdsum - Weirdsum2*Weirdsum3
		#print(C300Method2,C300Method3)
		
		#Method 4
		if nt == Nt-1:
			
			#Energy
			Earrayavg = np.average(Earray)
			E2arrayavg = np.average(E2array)
			EX0 = E2arrayavg - Earrayavg*Earrayavg
			
			
			
			#Magnetization
			Marrayavg = np.average(Marray)
			M2arrayavg = np.average(M2array)
			MX0 = M2arrayavg - Marrayavg*Marrayavg
			
			
			
			for j in range(nautotimes):
				#MWeirdsum = 0
				#MWeirdsum2 = 0
				#MWeirdsum3 = 0
				
				#EWeirdsum = 0
				#EWeirdsum2 = 0
				#EWeirdsum3 = 0
				
				#Vectorized sums
				MWeirdsum = np.sum(Marray[:nautotimes-j]*Marray[j:nautotimes])
				MWeirdsum2 = np.sum(Marray[:nautotimes-j])
				MWeirdsum3 = np.sum(Marray[j:nautotimes])
				
				EWeirdsum = np.sum(Earray[:nautotimes-j]*Earray[j:nautotimes])
				EWeirdsum2 = np.sum(Earray[:nautotimes-j])
				EWeirdsum3 = np.sum(Earray[j:nautotimes])
				
				#for i in range(nautotimes-j):
				#	MWeirdsum += Marray[i]*Marray[i+j]
				#	MWeirdsum2 += Marray[i]
				#	MWeirdsum3 += Marray[i+j]
					
				#	EWeirdsum += Earray[i]*Earray[i+j]
				#	EWeirdsum2 += Earray[i]
				#	EWeirdsum3 += Earray[i+j]
				#	
					
				MWeirdsum *= (1.0/(nautotimes-j))
				MWeirdsum2 *= (1.0/(nautotimes-j))
				MWeirdsum3 *= (1.0/(nautotimes-j))
				
				EWeirdsum *= (1.0/(nautotimes-j))
				EWeirdsum2 *= (1.0/(nautotimes-j))
				EWeirdsum3 *= (1.0/(nautotimes-j))
				#C300Method4 = (Weirdsum - Weirdsum2*Weirdsum3)/X0
				CMts[j] = (MWeirdsum - MWeirdsum2*MWeirdsum3)/MX0
				CEts[j] = (EWeirdsum - EWeirdsum2*EWeirdsum3)/EX0

			#print(C300Method4)
		

		
		
		#==========================================
		#What data do other get?
		#2-dimensional_ising_model
		#XT i ranges 0-100
		#Cv i range 0-0.6
		#ENergy pr link 0-(-1)
		
		#advancedlab1_Part1
		#E per site 0-(-2)
		#Cv specific 0-0.12
		#X 0-0.3
		#Nok kun E som er pr site?
		
		
		#AsherIsingModelReport
		#E pr site 0-(-2)
		#X 0-0.4
		#Cv 0-0.2
		
		
		#Chap12_Ising_Model_v04
		#E pr site 0-(-2)
		#Cv/N 0-4
		
		
		#Ising (1)
		#C/N 0-3
		#E/N 0-(-2)
		
		
		#Ising Model 1
		#E 0-(-1)
		#Cv 0-10
		#XT 0-300
		
		
		#ising_Model
		#E pr spin 0-(-4)
		#Cv 0-5000
		#XT pr spin 0-0.18
		
		
		#lecture6-stat_mech
		#Cv pr spin 0-2
		#
		
		#SCALAS2013_04_16-17
		#Cv 0-1.5
		
		#Student Ising Swarthmore
		#E fra 0-(-2)
		#Cv 0-0.2
		#X 0-0.1


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

####################################
#jeg kan prøve at se med 3x3 case og se om den calculate det samme som jeg gør på paper, og måske

if __name__ == "__main__":
	nx = 35
	ny = 35
	N = nx*ny
	#x = np.linspace(0,10,nx)
	#y = np.linspace(0,10,ny)
	T = 0.9
	dT = 0.1
	Nt = 40
	J = 1# J > 0 gives ferromagnetism
	ntest = 16 #Amount of Monte Carlo runs we do for EACH temperature. To weed out local minimum effects
	
	nautotimes = 2000
	
	#======================================================
	#Make (T,M) graph etc. For each temperature T, we add the final averaged out M,E,XT,Cv
	Ms = []
	Es = []
	XTs = []
	Cvs = []
	Ts = []
	
	#Autocorrelation
	CEts = np.zeros(nautotimes)
	CMts = np.zeros(nautotimes)
	
	
	#1 for save, 0 for no savefig
	SaveFig = 0
	
	MainScript(T,Ms,Es,XTs,Cvs,Ts,CEts,CMts)
	
	
	
	#=====================================================
	#Make plots
	fig2 = plt.figure(2)
	plt.plot(Ts,Ms)
	plt.xlabel("T")
	plt.ylabel("m")
	plt.title("Magnetization")
	fig2.savefig('Magnetization.png', bbox_inches='tight')
	plt.show()

	fig3 = plt.figure(3)
	plt.plot(Ts,Es)
	plt.xlabel("T")
	plt.ylabel("E")
	plt.title("Energy pr spin site")
	fig3.savefig('Energy.png', bbox_inches='tight')
	plt.show()

	fig4 = plt.figure(4)
	plt.plot(Ts,XTs)
	plt.xlabel("T")
	plt.ylabel("XT")
	plt.title("Isothermal Susceptibility")
	fig4.savefig('Susceptibility', bbox_inches='tight')
	plt.show()


	fig5 = plt.figure(5)
	plt.plot(Ts,Cvs)
	plt.xlabel("T")
	plt.ylabel("Cv")
	plt.title("Specific heat")
	fig5.savefig('Specificheat.png', bbox_inches='tight')
	plt.show()


	fig6 = plt.figure(6)
	tshere = [j for j in range(nautotimes)]
	plt.plot(tshere,CMts, color="red")
	plt.plot(tshere,CEts, color="black")
	plt.legend(["M","E"])
	plt.title("Autocorrelation")
	plt.xlabel("n timestep")
	plt.ylabel("C(n)")
	fig5.savefig('Autocorrelation.png', bbox_inches='tight')
	plt.show()

	#print(Ms)
	#print(Ts)
	#Lige nu der går den til M = +-625... 
	
	
	#===============================================================
	#Animation
	#anim = animation.FuncAnimation(fig, animate, frames=10, interval=100)
	# fig = plt.figure()
	# ax = fig.gca()
	# #plt.imshow(S, interpolation=None,vmin = -1, vmax = 1,cmap = "jet")
	# ax.set_title('Ising model, Metropolis algorithm')
	# fig.savefig('Ising.png', bbox_inches='tight')
	# #plt.colorbar()
	# plt.show()