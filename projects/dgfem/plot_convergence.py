import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#used if no GUI is available
matplotlib.use('Agg')
 
# exponential function x = 10^y
datax = [2, 8, 32, 128, 512]
datay = [0.14546, 0.040082, 0.0105124, 0.075133, 0.0078868]
datay2 = [0.16666667, 0.0416667, 0.0104167, 0.00260417, 0.000651042]
 
#convert x-axis to Logarithmic scale
plt.xscale("log")
plt.yscale("log")

plt.plot(datax,datay, label="discontinuous")
plt.plot(datax,datay2, label="classic")
plt.label()
plt.savefig("convergence1", bbox_inches='tight')