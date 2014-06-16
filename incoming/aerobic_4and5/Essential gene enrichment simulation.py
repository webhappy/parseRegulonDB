    ### simulation to create a range of possible enrichments of essential genes 
    ### if picked at random from the PEC database

import numpy as np
import matplotlib.pyplot as plt

# parameters
attempts = 10000 # Number of simulated attempts made
picked = 100 # Number of genes picked in each attempt
PEC_e = float(302) # Number of essential genes in PEC
PEC_ne = float(4190) # Number of non-essential genes in PEC
essential = PEC_e/(PEC_e+PEC_ne) # Proportion of essential genes in E. coli according to PEC
enrichments = [] # list of the proportions of essentials in each attempt
t8_correct = 87 # Number of essential genes (as per PEC) in my top 100 hits

# Simulation
for p in range(attempts):
    temp = 0
    for q in range(picked):
        if np.random.uniform() < essential:
            temp += 1
    enrichments.append(temp)
    
# Plotting
plt.figure()
num_bins = range(101)
n, bins, patches = plt.hist(enrichments, num_bins, facecolor='green')
plt.plot([t8_correct, t8_correct], [0, 1800], 'r-', linestyle='--', lw=2)
plt.savefig('/Users/Toro/Documents/Dropbox/Data/Cas9/Miseq/Analysis/PEC comparison/essential gene enrichment.png')
plt.show()