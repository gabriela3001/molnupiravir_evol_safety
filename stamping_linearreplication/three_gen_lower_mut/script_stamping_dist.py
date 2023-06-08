import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
from collections import Counter


m = int(sys.argv[1])

mrate_range = np.linspace(1e-6,1e-3,100)
mrate = mrate_range[m]

def three_generations(mrate):

    first_gen = np.random.binomial(2, mrate/3, 4)
    parents_second_gen = np.random.choice(first_gen, 20)
    second_gen = np.random.binomial(2, mrate/3, 20) + parents_second_gen
    parents_third_gen = np.random.choice(second_gen, 100)
    third_gen = np.random.binomial(2, mrate/3, 100) + parents_third_gen
    
    return(list(third_gen).count(1))
    
def stamping(mrate):

    first_gen = np.random.binomial(2, mrate, 100)
    
    return(list(first_gen).count(1))
    
three_generations_all = np.array([three_generations(mrate) for _ in range(5000000)]).flatten()
stamping_all = np.array([stamping(mrate) for _ in range(5000000)]).flatten()

output = '/home/labs/pilpel/gabril/COVID_merck/plos_revision/three_gen_lower_mut/'

with open(output + 'results_' + str(m) + '.txt', 'wb') as f:
    pickle.dump({'three_gen_result': Counter(three_generations_all), 'stamping_result': Counter(stamping_all)}, f)