from all_formulas import *
import numpy as np
import pickle
import numpy as np
import sys


mval = int(sys.argv[1])

mvals_heatmap = np.linspace(10000,25000,800)
a1vals_heatmap = np.linspace(7.7,10,800)
u_vals = np.logspace(-7,-1.5,10000)

T = 5
b = 7.61
a0 = 3
u = 1e-6
n = 1


all_results = []

for a1ind in range(len(a1vals_heatmap)):
    a1val = a1vals_heatmap[a1ind]

    index_u0 = np.searchsorted(u_vals, u)
    index_usafe = np.searchsorted(-totalY_wholeinfection(T, b, a0, a1val, 1-u, 1-u_vals[index_u0:], mvals_heatmap[mval], n),-totalY_wholeinfection(T, b, a0, a1val, 1-u, 1-u, mvals_heatmap[mval], n))

    all_results.append({'b':b, 'T':T, 'a0':a0, 'u0':u, 'm':mvals_heatmap[mval], 'n':n, 'a1':a1val,
                        'u_safe': u_vals[index_u0+index_usafe], 'u_ratio': u_vals[index_u0+index_usafe]/u_vals[index_u0]})
    

    
output = '/home/labs/pilpel/gabril/COVID_merck/final_runs/heatmap_usafe_ratio/treatment_peak_short/' 
    
with open(output + 'results_'+str(mval)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)