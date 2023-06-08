import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
from sklearn.model_selection import ParameterGrid
import pickle
import sys

def merck_euler(population,dt,b,q,m,n,a,beta):
    
    y1,y2,y3 = population
    
    dy      = np.empty(3)
    
    dy[0] = y1 + dt*(y1*(b*(q**(m+n))-a-beta*y3))
    dy[1] = y2 + dt*(y1*b*(q**m)*(1-(q**n)) + y2*(b*(q**m)-a-beta*y3))
    dy[2] = y3 + dt*((y1+y2)*(b*(1-q**m)) - a*y3)
    
    return(dy)
	
def simulation(init,b,a,ap,q,qd,m,n,d,tim,ttr):
    
    dt = 1e-3
    initial = [init[0],init[1], init[2]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q
    current_a = a

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b,current_q,m,n,current_a,d)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = ap
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = qd
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)
	
param = int(sys.argv[1])

with open('paramgrid_newparams_a03.txt', 'rb') as f:
    paramgrid_file = pickle.load(f)
    
paramgrid = paramgrid_file[param]

b = 7.61
q = 1-1e-6
m = paramgrid['m']
n = paramgrid['n']
a = 3
ap = paramgrid['ap']
#qd = paramgrid['qd']
tim = 5
#ttr = paramgrid['ttr']
beta = 1e-8
along_qd = list(1-np.logspace(-6,-5,20)) + [1-3e-6]
along_ttr = list(np.linspace(0,14,20)) + [5]


init_pop = [1,0,0]

all_results = []

for qd in along_qd:
    print(along_qd.index(qd))
    start_time = time.time()
    for ttr in along_ttr:
        print(along_ttr.index(ttr))
        control = simulation(init_pop,b,a,ap,q,q,m,n,beta,tim,ttr)
        with_treatment = simulation(init_pop,b,a,ap,q,qd,m,n,beta,tim,ttr)
        Y_control = np.sum(control[:,1])*1e-2
        Y_treatment = np.sum(with_treatment[:,1])*1e-2
        all_results.append({'m':m, 'n':n, 'beta': beta, 'ap':ap, 'qd':qd, 'ttr':ttr, 'Y_control':Y_control, 'Y_treatment':Y_treatment})

    print(time.time() - start_time)
	
output = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/lethal_defection_betaminus8/'
with open(output+'results_colormap_'+str(param)+'.txt', 'wb') as f:
    pickle.dump({'param':int(param), 'results':all_results},f)