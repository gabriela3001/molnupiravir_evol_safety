import numpy as np
import matplotlib.pyplot as plt
import time as tm
import sys
import pickle

param = int(sys.argv[1])
niter = int(sys.argv[2])

with open('paramgrid_minus_plus_strand.txt', 'rb') as f:
    param_grid = pickle.load(f)

mrate = param_grid[param]['mrate']
progression_synthesis = param_grid[param]['progression_synthesis']
ncells = param_grid[param]['ncells']


def one_cell(mrate, progression_synthesis):

    synthesis_minus_strand = np.random.binomial(2, mrate, progression_synthesis[1])
    synthesis_plus_strand = np.random.binomial(2, mrate, progression_synthesis[2])
    adding_up = synthesis_plus_strand + np.random.choice(synthesis_minus_strand, size = len(synthesis_plus_strand))
    
    return(list(adding_up))

def all_patients(npatients, ncells, mrate, progression_synthesis):
    all_p = []
    for p in range(npatients):
        by_patient = np.zeros(3)
        for n in range(ncells):
            res = one_cell(mrate, progression_synthesis)
            by_patient += np.array([res.count(0), res.count(1), res.count(2)])
        all_p.append(by_patient)
    return(all_p)

start_time = tm.time()
sim = all_patients(10000, ncells, mrate, progression_synthesis)
print(tm.time()-start_time)

output_dir = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/plus_minus_strand/plus_minus_simulation/'

with open('results_plus_minus_'+str(param) + '_' + str(niter)+ '.txt', 'wb') as f:
    pickle.dump(sim, f)