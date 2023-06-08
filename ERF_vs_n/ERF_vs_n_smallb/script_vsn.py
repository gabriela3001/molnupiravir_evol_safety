from all_formulas import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.model_selection import ParameterGrid
import sys

b = 6.61
q = 1-1e-6
a0 = 2
tim = 5
init_pop = [1,0]
qd = 1-3e-6

with open('paramgrid_vsn.txt', 'rb') as f:
    param_grid = pickle.load(f)
    
p = int(sys.argv[1])

mval = param_grid[p]['m']
a1val = param_grid[p]['a1']-1

n_values = np.linspace(1,29900-mval,50)

ERF_peak = []
ERF_exposure = []

for n in range(len(n_values)):
    simulation_data_notreatment = simulation([1,0],b,a0,a1val,q,q,mval,n_values[n],tim,5)
    simulation_data_treatment_peak = simulation([1,0],b,a0,a1val,q,qd,mval,n_values[n],tim,5)
    simulation_data_treatment_exposure = simulation([1,0],b,a0,a1val,qd,qd,mval,n_values[n],tim,0)
    
    ERF_peak.append(np.sum(simulation_data_treatment_peak[:,1])/np.sum(simulation_data_notreatment[:,1]))
    ERF_exposure.append(np.sum(simulation_data_treatment_exposure[:,1])/np.sum(simulation_data_notreatment[:,1]))
    
results = {'a1':a1val, 'm':mval, 'ERF_peak':ERF_peak, 'ERF_exposure':ERF_exposure}

output = '/home/labs/pilpel/gabril/COVID_merck/final_runs/ERF_vs_n_smallb/'

with open(output + 'results_' + str(p) + '.txt', 'wb') as f:
    pickle.dump(results, f)