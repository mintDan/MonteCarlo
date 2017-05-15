"""
Calc Monte Carlo pi


Dan Krog



"""

from matplotlib import pyplot
import numpy as np
import time
import timeit
import multiprocessing
from multiprocessing import Pool
import sys


#===============================================
#Timing decorator

def TimingCode(f):
	def timed(*args, **kw):
		ts = time.time()
		result = f()#(*args, **kw)
		te = time.time()
		print("Time difference based on decorator")
		print(f.__name__)
		print(te-ts)
		
		return result #Den return result fra RunMonteCarlo()
	
	#Denne her return en function, ikke et tal etc, men en function, important difference
	return timed #Den return result fra timed, som return result fra RunMonteCarlo()
	

	
@TimingCode
def RunMonteCarlo():
	"""
	Hello there!!!
	"""
	#inse = np.array([])
	inse = np.array([],dtype=bool)
	inside = 0
	
	#Måske er det muligt at bruge array[xg*xg+yg*yg <= rsqrd], i stedet for for loop
	
	
	for i in range(n):
		if xg[i]*xg[i]+yg[i]*yg[i] <= rsqrd:
			inside+=1
			inse = np.append(inse,i)
	return inse,inside

#Decorator fjerner pt docstring fra RunMonteCarlo() desværre...
#print(RunMonteCarlo.__doc__)
#http://stackoverflow.com/questions/1622943/timeit-versus-timing-decorator
#Noget med man kan lave en decorator, så den IKKE fjerner docstrings...

#print(RunMonteCarlo.__name__)
#Yea, og, den print at name er timed, det skal det jo ikke være, name skal være RunMonteCarlo!!!
#Så yea, ngoet med man skal bruge en rigtig wrapper... den tager sikkert __doc__ og __name__ etc og giver det videre correctly


#inse, inside = RunMonteCarlo()

#===========================================
#Manual timing, no decorator
#t0timeit = timeit.default_timer()
#t0 = time.clock()
#inse, inside = RunMonteCarlo()
#t1timeit = timeit.default_timer()
#t1 = time.clock()


# print("Time for Monte Carlo, seconds, using time module")
# print(t1-t0)
# print("Time for Monte Carlo, seconds, using timeit module with default_timer()")
# print(t1timeit-t0timeit)



#==========================================
#MP MC
#@TimingCode
def RunMonteCarloMP(): #normalt tager den n
	"""
	Kan bruges til apply, right? Fordi den ikke tager args
	Men ikke til pool.map
	"""
	
	n = int(10000000/4)
	count = 0
	for i in range(n):
		x = np.random.uniform()
		y = np.random.uniform()
		
        # if it is within the unit circle
		if x*x + y*y <= 1:
			count=count+1
	return count

	
	
def RunMonteCarloMP2(n):
	"""
	Hello there!!!
	"""
	#inse = np.array([])
	#inse = np.array([],dtype=bool)
	inside = 0
	
	xg = np.random.uniform(-1, 1, size=n)
	yg = np.random.uniform(-1, 1, size=n)
	
	b = 0.5
	r=np.sqrt(b**2)
	rsqrd = r*r
	
	#Måske er det muligt at bruge array[xg*xg+yg*yg <= rsqrd], i stedet for for loop
	
	
	
	for i in range(n):
		if xg[i]*xg[i]+yg[i]*yg[i] <= rsqrd:
			inside+=1

	return inside
	

def RunMonteCarloMP3():
	"""
	Hello there!!!
	"""
	#inse = np.array([])
	#inse = np.array([],dtype=bool)
	inside = 0
	
	#xg = np.random.uniform(-1, 1, size=n)
	#yg = np.random.uniform(-1, 1, size=n)
	
	#n = len(xg)
	
	b = 0.5
	r=np.sqrt(b**2)
	rsqrd = r*r
	
	#Måske er det muligt at bruge array[xg*xg+yg*yg <= rsqrd], i stedet for for loop
	
	
	
	for i in range(n):
		if xg[i]*xg[i]+yg[i]*yg[i] <= rsqrd:
			inside+=1

	return inside


	
def RunMonteCarloMP4(): #normalt tager den n
	"""
	Kan bruges til apply, right? Fordi den ikke tager args
	Men ikke til pool.map?
	ELler, map tager da også args, det er pool.map(f,iterable)
	"""
	
	
	n = int(10000000/4)
	
	try:
		print(lul)
	except:
		print("Couldnt print global variable")
	
	count = 0
	for i in range(n):
		x = np.random.uniform()
		y = np.random.uniform()
		
        # if it is within the unit circle
		if x*x + y*y <= 1:
			count+=1
	return count
	
	
	
	
	
#==========================================
#To plot the function
def PlotStuff():
	xr = np.linspace(a,b,n)
	def f1(x):
		return np.sqrt(r**2-xr**2)
	def f2(x):
		return -np.sqrt(r**2-xr**2)

	#pyplot.hold(True)
	pyplot.scatter(xg[::k],yg[::k],color='black')
	pyplot.scatter(xg[inse][::k],yg[inse][::k],color='black')
	pyplot.plot(xr,f1(xr),color='red', lw = 3)
	pyplot.plot(xr,f2(xr),color='red',lw = 3)
	#pyplot.axis([x0,x1,y0,y1])
	pyplot.xlabel('x')
	pyplot.ylabel('y')
	pyplot.title('Approximate $\pi$')
	pyplot.legend(['$\pi$ = {0:2f}'.format(Pi)])
	pyplot.savefig('MCpi.png', bbox_inches='tight')
	pyplot.show()

if __name__ == "__main__":
	
	#Grid
	x0 = -1
	x1 = 1
	y0 = -1
	y1 = 1
	#Circle
	a = -0.5
	b = 0.5
	r=np.sqrt(b**2)
	rsqrd = r*r

	k = 100
	n = k*1000

	#x = np.linspace(x0,x1,n)
		
	#==============================================
	#Make guesses/throw points
	xg = np.random.uniform(x0, x1, size=n)
	yg = np.random.uniform(y0, y1, size=n)
	
	
	print("Size of xg with sys.getsizeof(), giving bytes")
	print(sys.getsizeof(xg))
	print("Size of xg with .__sizeof__(), giving bytes")
	print(xg.__sizeof__())
	print("Size pr element with sys.getsizeof()")
	print(sys.getsizeof(xg)/n)
	#===============================================
	#Run MC script
	
	inse, inside = RunMonteCarlo()
	
	
	#===================================================
	#calculating probabilities and pi$		
	P = inside/n
	print('Probability to land inside circle')
	print(P)


	print('Analytical area')
	print(np.pi*r**2)


	# Pi = P*(b-a)**2/(r**2)
	# print('Approximated Pi')
	# print('Pi={}'.format(Pi))

	Pi = P*(x1-x0)*(y1-y0)/(r**2)
	print('Approximated Pi')
	print('Pi={}'.format(Pi))
	
	
	#=====================================================
	#To see if they spot a global
	lul = "Hi from Global!"
	
	#=====================================================
	#Plot 
	#PlotStuff()
	numberp = multiprocessing.cpu_count()
	print("Number of CPUs: {0:1d}".format(numberp))
	# Nummber of points to use for the Pi estimation
	n = 10000000
	# iterable with a list of points to generate in each worker
	# each worker process gets n/np number of points to calculate Pi from
	split_n =[int(n/numberp) for i in range(numberp)]
	#Create the worker pool
    # http://docs.python.org/library/multiprocessing.html#module-multiprocessing.pool
	pool = Pool(processes=numberp)   
	# parallel map
	#count=pool.map(RunMonteCarloMP, part_count)
	count = pool.apply_async(RunMonteCarloMP)
	
	#print("Esitmated value of Pi:: ", sum(count)/(n*1.0)*4)
	ntemp = int(10000000/4)
	print(count.get()/(ntemp)*4)
	
	
	#Med apply_async skal man også bruge .get() bagefter
	#inside = pool.apply_async(RunMonteCarloMP2)
	
	#====================
	#RunMonteCarloMP2 kan ikke se GlobalTest når jeg kører den fra pool.map()
	#Så, måske skal man ikke have global arrays med Multiprocessing, eller noget andet...
	#
	#GlobalTest = 23
	
	#===================================================================================================
	#Okay, så function i pool.map() kan ikke se guesses... Så, xg,yg skal defines inde i function i guess
	#Så, med pool.map kan den ikke se global varialbes eller hvad? KUn ting der bliver passed til function, og defined inde i function?
	#Ah, inside er jo nok en list her... en list med 4 numbers, fordi arg er part_count som har 4 items
	
	#man kunne også lave split_n = [(n1,n2),(n2,n3),(n3,n4)]
	#Så kan man nemlig lavet for loop fra n1,n2 inde i funktion...
	#Men, så nytter det stadig ikke noget at man ikke kan se global arrays...
	#Ahh... 
	#Men, btw, kan den return numpy array, RunMonteCarloMP2?
	#Bare any numpy array?
	#Fordi jeg kunne jo også lave, pool.map(RunMonteCarloMP2,[u[n1:n2],u[n2:n3],u[n3:n4])
	#Så jeg GIVER de arrays til RunMonteCarloMP2 som den skal regne på...
	#Også kan man jo bare tage n=len(array) inde i function
	#MEn damn, der er to arrays...
	#Mon at (xg,yg) den copy de to arrays, aka tager lang tid mon?
	#Jeg skal pass TO arrays til RunMonteCarloMP3
	
	#Så altså, det kan være at apply_async kan se global variables??
	#map kan vidst ikke
	
	inside = pool.map(RunMonteCarloMP2, split_n)
	print(inside)
	inside = np.sum(inside)
	
	
	P = inside/n
	print('Probability to land inside circle')
	print(P)
	
	Pi = P*(x1-x0)*(y1-y0)/(r**2)
	print('Approximated Pi')
	print('Pi={}'.format(Pi))
	
	
	#===========================================================
	#Use RunMonteCarloMP4()
	#Tager ingen args
	#pool.apply_async
	#Fungerer IKKE med global variables, den kan IKKE finde global variables.
	#Heller ikke med global lul inde i function, kan den finde lul 
	#heller ikke n
	
	#numberp = multiprocessing.cpu_count()
	print("Number of CPUs: {0:1d}".format(numberp))
	# Nummber of points to use for the Pi estimation
	ntotal = 10000000
	# iterable with a list of points to generate in each worker
	# each worker process gets n/np number of points to calculate Pi from

	n = int(ntotal/numberp)
	#Create the worker pool
    # http://docs.python.org/library/multiprocessing.html#module-multiprocessing.pool
	
	
	
	pool1 = Pool(processes=numberp)   
	# parallel map
	#count=pool.map(RunMonteCarloMP, part_count)

	count = [pool1.apply_async(RunMonteCarloMP4) for _ in range(numberp)]
	Counting = [res.get() for res in count]
	print("Result from each pool/worker")
	print(Counting)
	SumIt = np.sum(Counting)
	
	#print("Esitmated value of Pi:: ", sum(count)/(n*1.0)*4)
	print(SumIt/ntotal*4)
	
	
	
	