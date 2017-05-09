from matplotlib import pyplot
import numpy as np
import time
import timeit
import multiprocessing
from multiprocessing import Pool


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
	Hello there!!!
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
	#Plot 
	#PlotStuff()
	np = multiprocessing.cpu_count()
	print("Number of CPUs: {0:1d}".format(np))
	# Nummber of points to use for the Pi estimation
	n = 10000000
	# iterable with a list of points to generate in each worker
	# each worker process gets n/np number of points to calculate Pi from
	part_count=[int(n/np) for i in range(np)]
	#Create the worker pool
    # http://docs.python.org/library/multiprocessing.html#module-multiprocessing.pool
	pool = Pool(processes=np)   
	# parallel map
	#count=pool.map(RunMonteCarloMP, part_count)
	count = pool.apply_async(RunMonteCarloMP)
	
	#print("Esitmated value of Pi:: ", sum(count)/(n*1.0)*4)
	ntemp = int(10000000/4)
	print(count.get()/(ntemp)*4)