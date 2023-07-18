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
executables = [3, 4, 5, 6, 7, 8, 9, 10]

#domain info
h = [0.707596, 0.565214, 0.388975, 0.288115, 0.197503, 0.136591, 0.104136, 0.0703841, 0.0489181, 0.0364704]
n_cells = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
h_s = [i**2 for i in h]
h_2 = [i/12 for i in h_s]

    
def extract(run_name, measurements_dir, n_executables):
    measurements_dict = {}
    n_cells_list = []
    c_invs = []
    c_sigmas = []

    for exec in executables:

        run_and_exec = f'{run_name}_{exec}'

        for filename in os.listdir(run_and_exec):
            f = os.path.join(run_and_exec, filename)

            
            #accuire info about measurement
            filename_split = filename.split("_")

            c_inv = float(filename_split[1])
            c_sigma = float(filename_split[2])
            n_cells = int(filename_split[0])
            error_kind = filename_split[3][:-4]

            if c_sigma not in c_sigmas:
                c_sigmas.append(c_sigma)
            if c_inv not in c_invs:
                c_invs.append(c_inv)
            if n_cells not in n_cells_list:
                n_cells_list.append(n_cells)
            

            error = 100.0

            #read error from file
            with open(f) as myfile:
                first_line = myfile.readline()
                error = float(first_line)

            new_tuple = (int(exec), float(c_inv), float(c_sigma), int(n_cells), error_kind)
            measurements_dict[new_tuple] = error

        c_invs.sort()
        c_sigmas.sort()
        n_cells_list.sort()

    return measurements_dict, c_invs, c_sigmas, n_cells_list

#plot L2 error against O(h)
def plot_convergence_O_h(c_inv, c_sigma, n_cells_list, error, exec1, label1, exec2 = 0, label2 = ""):

    x = n_cells_list
    y = []
    y2 = [] 

    for n_cells in n_cells_list:
        l2_error = measurements_dict[(exec1, c_inv, c_sigma, n_cells, error)]
        y.append(l2_error)
    
    if (exec2 != 0):
        for n_cells in n_cells_list:
            l2_error = measurements_dict[(exec2, c_inv, c_sigma, n_cells, error)]
            y2.append(l2_error)
    
    X = np.array(x)
    Y = np.array(y)

    #convert x-axis to Logarithmic scale
    plt.figure(figsize=(15,10))
    fig = plt.figure(dpi=300)
    plt.xscale("log")
    plt.yscale("log")

    # Label the axes
    plt.xlabel('n [# cells in mesh]')
    plt.ylabel('L2 error of solution')
    plt.title(r'Convergence of manifactured solution $\mathbf{1 + x + y^2} $')

    plt.scatter(x, y, label=label1, marker = "p", color ="darkred")
    if (exec2 != 0):
        plt.scatter(x, y2, label=label2, marker = "p", color ="darkblue")
        
    plt.plot(x, h_2, label=r'$\mathcal{O}(h^2)$', linestyle="dotted", color ="green")

    plt.legend()
    plt.savefig(f'plots/{run_name}/L2_{c_inv}_{c_sigma}_{exec}_h.png')


measurements_dict, c_invs, c_sigmas, n_cells_list = extract(run_name, measurements_dir)


measurements_dict, c_invs, c_sigmas, n_cells_list = extract(run_name, measurements_dir)

print("C_invs:")
print(c_invs)
print("C_sigmas")
print(c_sigmas)
print("N_cells:")
print(n_cells_list)






 
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