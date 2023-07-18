import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os
import sys
import math

#used if no GUI is available
matplotlib.use('Agg')
 
# exponential function x = 10^y

h = [0.707596, 0.565214, 0.388975, 0.288115, 0.197503, 0.136591, 0.104136, 0.0703841, 0.0489181, 0.0364704]

datax = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
datay = [0.0295876, 0.0197537, 0.00946608, 0.0063399, 0.00341458, 0.00218018, 0.00121912, 0.000661853, 0.000386771, 0.00019453]

h_s = [i**2 for i in h]
h_2 = [i/12 for i in h_s]

params = np.polyfit(datax, h_2, 1)
predicter = np.poly1d(params)
print(params)

h_2_poly = [predicter(i) for i in datax]

data_ref = []
ref0 = 0.029
for i in range(len(datay)):
    new_num = ref0 / (2**i)
    data_ref.append(new_num)
    

datay2 = [0.16666667, 0.0416667, 0.0104167, 0.00260417, 0.000651042]
 
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