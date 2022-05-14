#this file is used to visualize legendre polynomials in 2D on tht unit square

from sys import argv

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from scipy.special import legendre

#used if no GUI is available
matplotlib.use('Agg')

fig = plt.figure(figsize=(30, 30))

x = np.arange(-1,1,0.1)
y = np.arange(-1,1,0.1)

X, Y = np.meshgrid(x,y)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


for i in range(3):
    for j in range(3):

        legendre_x = legendre(i)
        legendre_y = legendre(j)

        z = legendre_x(X) * legendre_y(Y)

        ax = fig.add_subplot(3,3,i*3+j+1, projection='3d')
        ax.plot_surface(X, Y, z, cmap=cm.coolwarm)
        ax.set_title(r'Legendre polynomial L_{i}(x) * L_{j}(y)')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

plt.savefig("legendre", bbox_inches='tight')
