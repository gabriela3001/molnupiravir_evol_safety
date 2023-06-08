import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time as tm
import pickle
import sys

with open('paramgrid_ERF_epistasis.txt', 'rb') as f:
    paramgrid = pickle.load(f)

#

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
    
    dt = 1e-4
    initial = [init[0],init[1],init[2],init[3]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q0
    current_a = a0

    for ngen in np.arange(150/dt):
        if ngen % 100 == 0:
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


all_landscapes = [[7.61,7.61,7.61,7.61],
                  [7.61,7.61,7.61,7.61*1.01],
                  [7.61,7.61*1.01,7.61,7.61],
                  [7.61,7.61*1.01,7.61*1.01,7.61],
                  [7.61,7.61*1.01,7.61*1.01,7.61*1.01],
                  [7.61,7.61*0.97,7.61*0.99,7.61*1.01]]




param = int(sys.argv[1])

b = 7.61
q0 = 1-5e-6
m = paramgrid[param]['m']
a0 = 3
a1 = paramgrid[param]['ap']
tim = 5
along_qd = list(1-np.logspace(-6,-5,20)*5) + [1-3e-6*5]
along_ttr = list(np.linspace(0,14,20)) + [5]
b_array = all_landscapes[paramgrid[param]['landscape_index']]


init = [1,0,0,0]

all_results = []

start_time = tm.time()

for q1 in along_qd:
    print(along_qd.index(q1))
    for ttr in along_ttr:
        print(along_ttr.index(ttr))
        control = simulation(init,b_array,a0,a1,q0,q0,m,tim,ttr)*1e-3
        with_treatment = simulation(init,b_array,a0,a1,q0,q1,m,tim,ttr)*1e-3        
        all_results.append({'m':m, 'ap':a1, 'qd':q1, 'ttr':ttr, 'Y_control':np.sum(control, axis = 0), 
                            'Y_treatment':np.sum(with_treatment, axis = 0)})
        
end_time = tm.time()-start_time

output = '/home/labs/pilpel/gabril/COVID_merck/review/higher_u0/DM_u0_5_10e6/'

with open(output + 'results_'+str(param)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)