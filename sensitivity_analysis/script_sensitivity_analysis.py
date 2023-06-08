import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
from sklearn.model_selection import ParameterGrid
import pickle
import sys

def merck_euler(population,dt,b,q,m,n,a):

    mutation = 1
    
    y1,y2 = population
    
    dy      = np.empty(2)

    dy[0] = y1 + dt*(b*y1*(q**(m)*(q**n))-a*y1)
    dy[1] = y2 + dt*(mutation*(b*y1*(q**m)*(1-q**n)+b*y2*(q**m)-a*y2))
    
    return(dy)

def simulation(init,b,a,ap,q,qd,m,n,tim,ttr):
    
    dt = 1e-3
    initial = [init[0],init[1]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q
    current_a = a

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b,current_q,m,n,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = ap
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = qd
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)
    
with open('paramgrid_sensitivity_analysis.txt', 'rb') as f:
    paramgrid_file = pickle.load(f)
with open('paramgrid_b_a_values.txt', 'rb') as f:
    b_values, a_values = pickle.load(f)

param = int(sys.argv[1])  
paramgrid = paramgrid_file[param]

b = b_values[paramgrid['param']]
q = 1-1e-6
m = paramgrid['m']
n = 87
a = 3
ap = a_values[paramgrid['param']][paramgrid['aind']]
tim = paramgrid['param'][1]
ttr = 3
qd = 1-3e-6

init_pop = [1,0]


control = simulation(init_pop,b,a,ap,q,q,m,n,tim,ttr)
with_treatment = simulation(init_pop,b,a,ap,q,qd,m,n,tim,ttr)
Y_control = np.sum(control[:,1])*1e-2
Y_treatment = np.sum(with_treatment[:,1])*1e-2
all_results = {'m':m, 'n':n, 'ap':ap, 'qd':qd, 'ttr':ttr, 'Y_control':Y_control, 'Y_treatment':Y_treatment, 'b': b, 'ap':ap, 'tim' = tim}



output = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/sensitivity_analysis/'
with open(output+'results_othervariants_'+str(param)+'.txt', 'wb') as f:
    pickle.dump({'param':int(param), 'results':all_results},f)