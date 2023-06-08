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

def simulation(init,b,a0,a1,q0,q1,m,n,tim,ttr,ttr_end):
    
    dt = 1e-3
    initial = [init[0],init[1]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q0
    current_a = a0

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b,current_q,m,n,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = a1
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = q1
            
        # Merck treatment stops
        if ngen > ttr_end/dt:
            current_q = q0
            
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
along_qd = list(1-np.logspace(-6,-5,20)) + [1-3e-6]
along_ttr = list(np.linspace(0,14,20)) + [5]

init_pop = [1,0]

all_results = []

for qd in along_qd:
    print(along_qd.index(qd))
    for ttr in along_ttr:
        print(along_ttr.index(ttr))
        control = simulation(init_pop,b,a,ap,q,q,m,n,tim,ttr,ttr+5)
        with_treatment = simulation(init_pop,b,a,ap,q,qd,m,n,tim,ttr,ttr+5)
        Y_control = np.sum(control[:,1])*1e-2
        Y_treatment = np.sum(with_treatment[:,1])*1e-2
        
        Y_control_after = np.sum(control[:,1][int((ttr+5)*1e3):])*1e-2
        Y_treatment_after = np.sum(with_treatment[:,1][int((ttr+5)*1e3):])*1e-2
        
        all_results.append({'m':m, 'n':n, 'ap':ap, 'qd':qd, 'ttr':ttr, 'Y_control':Y_control, 'Y_treatment':Y_treatment,
                           'Y_control_after':Y_control_after, 'Y_treatment_after':Y_treatment_after})
        

output = '/home/labs/pilpel/gabril/COVID_merck/review/treatment_length/'
with open(output+'results_colormap_'+str(param)+'.txt', 'wb') as f:
    pickle.dump({'param':int(param), 'results':all_results},f)
