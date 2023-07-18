import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os
import sys
import math

#used if no GUI is available
matplotlib.use('Agg')
 
#directories info
run_name = sys.argv[1]
measurements_dir="measurements"
plots_dir="plots"

#domain info
h = [0.707596, 0.565214, 0.388975, 0.288115, 0.197503, 0.136591, 0.104136, 0.0703841, 0.0489181, 0.0364704]
n_cells = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
h_s = [i**2 for i in h]
h_2 = [i/12 for i in h_s]

    

 
#convert x-axis to Logarithmic scale
plt.figure(figsize=(15,10))
fig = plt.figure(dpi=300)
plt.xscale("log")
plt.yscale("log")

# Label the axes
plt.xlabel('n [# cells in mesh]')
plt.ylabel('L2 error of solution')
plt.title(r'Convergence of manifactured solution $\mathbf{1 + x + y^2} $')

plt.plot(datax, data_ref, label=r'$\mathcal{O}(n^{-2})$', linestyle="dotted", color ="black")
plt.plot(datax, h, label="h", linestyle="dotted", color ="grey")
plt.scatter(datax, datay, label="discontinuous", marker = "p", color ="darkred")
plt.plot(datax, h_2, label=r'$\mathcal{O}(h^2)$', linestyle="dotted", color ="green")

plt.legend()
plt.savefig("convergence_16_juli_1", bbox_inches='tight')