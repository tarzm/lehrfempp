from sys import argv

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#used if no GUI is available
matplotlib.use('Agg')

fig = plt.figure(figsize=(30, 30))

if len(argv) < 2:
    print('usage: python polytopic_plot.py mesh_function.txt -s')
    print('options: -s is simple mesh_function plot')
    print('simple is default')
    exit(-1)

filename = argv[1]

input = np.loadtxt(filename)

x = input[:,0]
y = input[:,1]
z = input[:,2]


X, Y = np.meshgrid(x,y)

png_name = filename[:-3]

ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_trisurf(x, y, z, cmap=cm.coolwarm)
ax.set_title(r'Global mesh function')
ax.set_xlabel("X")
ax.set_ylabel("Y")

plt.savefig(png_name, bbox_inches='tight')