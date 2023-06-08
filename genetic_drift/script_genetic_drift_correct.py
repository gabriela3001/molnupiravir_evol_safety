import numpy as np
import matplotlib.pyplot as plt
import time as tm
import numpy as np
import sys
import pickle

p = int(sys.argv[1])

with open('paramgrid_genetic_drift_short.txt', 'rb') as f:
	param_grid = pickle.load(f)

update_vectors = [[1,0], # x replicates without a mutation
                  [0,1], # y replicated without a mutation
                  [-1,0], # x replicates with lethal mutation
                  [0,1], # x replicates with VoC mutation
                  [-1,0], # x cleared by the immune system
                  [0,-1], # y replicates with lethal mutation
                  [0,-1] # y cleared by the immune system
                 ]
cat_names = ['x', 'y']		

def calculate_reactions(population, b, m, n, a, q):
    
    x, y = population
    
    rates = [x*b*(q**(m+n)),
             y*b*(q**m),
             x*b*(1-q**m),
             x*b*(1-q**n),
             a*x,
             y*b*(1-q**m),
             a*y]
    
    return(np.array(rates))
	

## initialize values

paramgrid = param_grid[p]

b = 7.61
q0 = 1-1e-6
m = paramgrid['m']
n = 87
a0 = 3
a1 = paramgrid['ap']
qd = paramgrid['qd']
tim = 5
ttr = 5
#along_qd = list(1-np.logspace(-6,-5,20)) + [1-3e-6]
#along_ttr = list(np.linspace(0,14,20)) + [5]

all_results = []
all_cumul = []

for nrun in range(1000):

    print(nrun)
    pop_int = np.zeros(len(cat_names))
    pop_int[0] = 5
    pop_int[1] = 0
    
    population = pop_int.copy()

    start_time = tm.time()
    
    times = [0]
    
    reaction_rates = calculate_reactions(population, b, m, n, a0, q0)
    
    counter = 0
    time = 0
    
    appearance_MT = False
    
    
    cumul_load = np.zeros(2)
    
    while time < 15:

        counter += 1

        # Reaction Time        
        tau = 1e-4
        time += tau
        
        # what happened
        events_happened = np.array(np.random.poisson([tau*reaction_rates[i] for i in range(len(reaction_rates))]))
        update_happened = np.array([np.array(update_vectors[i])*events_happened[i] for i in range(len(reaction_rates))])
        
        population = np.array(population) + np.sum(update_happened, axis = 0)
        population[population < 0] = 0
            

        # Update Reactions Rates
        if time < tim:
            a = a0
        else:
            a = a1
            
        ## update mutation rates
        if time < ttr:
            q = q0
        else:
            q = qd
            
        reaction_rates = calculate_reactions(population, b, m, n, a, q)
        
        cumul_load += population*tau
        
    end_time = tm.time() - start_time

    print(end_time)
    
    all_cumul.append(cumul_load)
    
    
result_dict = {'m':m, 'a1':a1, 'q1': qd, 'result':all_cumul}
    
output = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/genetic_drift_onepanel/'

with open('results_genetic_drift_'+str(p) + '.txt', 'wb') as f:
    pickle.dump(result_dict, f)