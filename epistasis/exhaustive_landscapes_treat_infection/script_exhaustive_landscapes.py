import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time as tm
import pickle
import sys

def merck_euler(population,dt,b_array,q,m,a):
    
    b00,b01,b10,b11 = b_array

    y00, y01, y10, y11 = population
    
    B00 = b00*(q**m)
    B01 = b01*(q**m)
    B10 = b10*(q**m)
    B11 = b11*(q**m)
    
    dy      = np.empty(4)

    dy[0] = y00 + dt*((B00*(q**2)*y00)
                     +(B01*q*(1-q)*y01)
                     +(B10*q*(1-q)*y10)
                     +(B11*((1-q)**2)*y11)
                     +(-a*y00))
    
    dy[1] = y01 + dt*((B00*q*(1-q)*y00)
                     +(B01*(q**2)*y01)
                     +(B10*((1-q)**2)*y10)
                     +(B11*q*(1-q)*y11)
                     +(-a*y01))
    
    dy[2] = y10 + dt*((B00*q*(1-q)*y00)
                     +(B01*((1-q)**2)*y01)
                     +(B10*(q**2)*y10)
                     +(B11*q*(1-q)*y11)
                     +(-a*y10))
    
    dy[3] = y11 + dt*((B00*((1-q)**2)*y00)
                     +(B01*q*(1-q)*y01)
                     +(B10*q*(1-q)*y10)
                     +(B11*(q**2)*y11)
                     +(-a*y11))

    return(dy)

def simulation(init,b_array,a0,a1,q0,q1,m,tim,ttr):
    
    dt = 1e-3
    initial = [init[0],init[1],init[2],init[3]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q0
    current_a = a0

    for ngen in np.arange(150/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b_array,current_q,m,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = a1
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = q1
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)
    
### preparing landscape
adv = 1.01
disadv = 0.99

def create_barray_landscape(b, landscape):
    barray = np.zeros(4)
    barray[0] = b
    
    for i in range(3):
        if landscape[i] == -1:
            barray[i+1] = disadv*b
        elif landscape[i] == 0:
            barray[i+1] = b
        elif landscape[i] == 1:
            barray[i+1] = adv*b
            
    return(barray)
    
with open('paramgrid_exhaustive_landscapes.txt', 'rb') as f:
  paramgrid = pickle.load(f)
    
param = int(sys.argv[1])

b = 7.61
q0 = 1-1e-6
m = paramgrid[param]['m']
a0 = 3
a1 = paramgrid[param]['a1']
tim = 5
q1_range = 1-np.logspace(-6,-3,100)
ttr = 0
landscape = paramgrid[param]['landscape']
#landscape = [0,1,1]
b_array = create_barray_landscape(b, landscape)
print(b_array)

initial = [1,0,0,0]

results = []

start_time = tm.time()

for q1 in q1_range:
    

    with_treatment = simulation(initial,b_array,a0,a1,q0,q1,m,tim,ttr)
    
    results.append({'landscape':landscape, 'a1':a1, 'm':m, 'u1':1-q1, 'Y00':np.sum(with_treatment[:,0])*1e-2, 
                      'Y01':np.sum(with_treatment[:,1])*1e-2, 'Y10':np.sum(with_treatment[:,2])*1e-2,
                      'Y11':np.sum(with_treatment[:,3])*1e-2})
                      
print(tm.time() - start_time)
                      
                      
output = '/home/labs/pilpel/gabril/COVID_merck/review/epistasis/exhaustive_landscapes_treat_infection/'

with open(output + 'results_'+str(param)+'.txt', 'wb') as f:
    pickle.dump(results, f)