import numpy as np
import sys
import pickle

def merck_euler_differentb(population,dt,b1,b2,q,m,n,a):
    mutation = 1
    x,y = population
    dy      = np.empty(2)
    dy[0] = x + dt*(b1*x*(q**(m)*(q**n))-a*x)
    dy[1] = y + dt*(mutation*(b1*x*(q**m)*(1-q**n)+b2*y*(q**m)-a*y))
    return(dy)
    
def simulation_precise_differentb(init,b1,b2,a,ap,q,qd,m,n,tim,ttr):
    
    dt = 1e-3
    initial = [init[0],init[1]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q
    current_a = a

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler_differentb(population,dt,b1,b2,current_q,m,n,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = ap
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = qd
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)
    
param = int(sys.argv[1])

with open('paramgrid_otheradvantages_halfperc.txt', 'rb') as f:
    param_file = pickle.load(f)
    
paramgrid = param_file[param]
    
b1 = 7.61
b2 = 7.61*paramgrid['adv']
a = 3
a1 = paramgrid['a1']
q = 1-1e-6
tim = 5
ttr= paramgrid['ttr']
m = paramgrid['m']
n=1
qd_range = 1-np.logspace(-6,-3,100)
init = [1, 0]


all_results = []
for qd in qd_range:
    #print(nval)
    sim = np.sum(simulation_precise_differentb(init,b1,b2,a,a1,q,qd,m,n,tim,ttr)[:,1])*1e-2
    all_results.append(sim)
    
output = '/home/labs/pilpel/gabril/COVID_merck/review/other_advantages_halfperc/'

with open(output+'results_advb_'+str(param)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)
