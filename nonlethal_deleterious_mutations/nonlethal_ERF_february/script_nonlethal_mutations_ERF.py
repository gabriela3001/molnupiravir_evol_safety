import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
from sklearn.model_selection import ParameterGrid
import pickle
import sys

def merck_euler(population,dt,b,b1,q,m,n1,n2,n3,a):
    
    x, y1,y2,y3 = population
    
    dy      = np.empty(4)
    
    M = m+n1+n2+n3
    
    dy[0] = x + dt*(x*(b*q**(M)-a))
    dy[1] = y1 + dt*(x*b*q**(M-n1)*(1-q**n1) + y1*(b*q**(M-n1)-a))
    dy[2] = y2 + dt*(x*b*q**(M-n2)*(1-q**n2) + y2*(b1*q**(M-n2)-a))
    dy[3] = y3 + dt*(x*b*(q**m)*(1-q**n3)+y1*b*(q**m)*(1-q**(n2+n3))+y2*b1*(q**m)*(1-q**(n1+n3)) + y3*(b1*(q**m)-a))
    
    return(dy)
    
def simulation(init,b,b1,a,ap,q,qd,m,n1,n2,n3,tim,ttr):
    
    dt = 1e-4
    initial = [init[0],init[1], init[2], init[3]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q
    current_a = a

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b,b1,current_q,m,n1,n2,n3,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = ap
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = qd
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)
       

with open('paramgrid_nonlethalmutations.txt', 'rb') as f:
    param_grid_list = pickle.load(f)
    
param = int(sys.argv[1])

param_grid = param_grid_list[param]

b = 7.61
q = 1-1e-6
m = param_grid['m']
n1 = 87
a = 3

if not param_grid['n2']:
  n2 = 0
  n3 = 0
else:
  n2 = 29900-87-100-3000-m
  n3 = 100


ap = param_grid['a1']
b1 = param_grid['b1']*b
tim = 5
ttr = 2

qd = 1-3e-6

init = [1,0,0,0]

all_results = []

control = simulation(init,b,b1,a,ap,q,q,m,n1,n2,n3,tim,ttr)
with_treatment = simulation(init,b,b1,a,ap,q,qd,m,n1,n2,n3,tim,ttr)
Y_control = np.sum(control[:,1])*1e-2 + np.sum(control[:,3])*1e-2
Y_treatment = np.sum(with_treatment[:,1])*1e-2 + np.sum(with_treatment[:,3])*1e-2
all_results.append({'m':m, 'n1':n1, 'ap':ap, 'qd':qd,'b1': b1,'n2':n2, 'n3':n3, 'ttr':ttr, 'Y_control':Y_control, 'Y_treatment':Y_treatment})
     
output = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/nonlethal_ERF_february/'
with open(output+'results_colormap_'+str(param)+'.txt', 'wb') as f:
    pickle.dump({'param':int(param), 'results':all_results},f)