
import numpy as np
import pickle
import numpy as np
import sys


mval = int(sys.argv[1])

a1vals_heatmap = np.linspace(7.7,10,800)
mvals_heatmap = np.linspace(10000,25000,800)
u_vals = np.logspace(-7,-3.5,10000)

T = 5
b = 7.61
a0 = 3
u = 1e-6
n = 1

## analytical functions for treatment at exposure

def plot_V(b, a0, a1, m, n, q, T):
    n1 = (a1-a0)*np.exp(T*(b*q**m-a0))
    d1 = (b*q**m-a0)*(a1-b*q**m)
    return(n1/d1)

def plot_X(b, a0, a1, m, n, q, T):
    n1 = (a1-a0)*np.exp(T*(b*q**(m+n)-a0))
    d1 = (b*q**(m+n)-a0)*(a1-b*q**(m+n))
    return(n1/d1)   

def plot_Y(b, a0, a1, m, n, q, T):
    X = plot_X(b, a0, a1, m, n, q, T)
    V = plot_V(b, a0, a1, m, n, q, T)
    return(V-X)



all_results = []

for a1ind in range(len(a1vals_heatmap)):
    print(a1ind)
    a1val = a1vals_heatmap[a1ind]

    index_u0 = np.searchsorted(u_vals, 1e-6)
    index_usafe = np.searchsorted(-plot_Y(b, a0, a1val, mvals_heatmap[mval], 1, 1-u_vals[index_u0:], T),-plot_Y(b, a0, a1val, mvals_heatmap[mval], 1, 1-1e-6, T))

    all_results.append({'b':b, 'T':T, 'a0':a0, 'u0':u, 'm':mvals_heatmap[mval], 'n':n, 'a1':a1val,
                        'u_safe': u_vals[index_u0+index_usafe], 'u_ratio': u_vals[index_u0+index_usafe]/u_vals[index_u0]})
    
output = '/home/labs/pilpel/gabril/COVID_merck/final_runs/heatmap_usafe_ratio/treatment_exposure_short/' 
    
with open(output + 'results_'+str(mval)+'.txt', 'wb') as f:
    pickle.dump(all_results, f)
                       