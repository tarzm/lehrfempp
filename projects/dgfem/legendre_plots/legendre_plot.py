#this file is used to visualize legendre polynomials in 2D on tht unit square

from sys import argv

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
import numpy as np
plt.rcParams.update(plt.rcParamsDefault)
from string import Template


from scipy.special import legendre

#used if no GUI is available
matplotlib.use('Agg')


fig = plt.figure(figsize=(30, 30))

x = np.linspace(-1,1,1000)
y = np.linspace(-1,1,1000)

X, Y = np.meshgrid(x,y)

# plt.rc('text', usetex=True)
# plt.rc('font', family='Latin Modern')

mathbf = '\mathbf'


for i in range(2):
    for j in range(2):

        legendre_x = legendre(i)
        legendre_y = legendre(j)

        z = legendre_x(X) * legendre_y(Y)



        ax = fig.add_subplot(3,3,i*3+j+1, projection='3d')
        ax.plot_surface(X, Y, z, cmap=cm.coolwarm)
        title = Template(r'$dollar\mathbf{\hat{\textit{L}}}_$iidx$dollar').substitute(iidx=i, dollar="$")
        ax.set_title(f'$\mathit{{\hat{{\Phi}}}}_{{{i*3+j}}}(\mathbf{{\hat{{x}}}}) = \mathbf{{\mathit{{\hat{{L}}}}}}_{{{i}}}(\hat{{x}}) \mathbf{{\mathit{{\hat{{L}}}}}}_{{{j}}}(\hat{{y}})$', fontsize=20)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

plt.savefig("legendre_new_1", bbox_inches='tight')
