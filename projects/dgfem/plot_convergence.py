import numpy as np
import os
import sys

#directories info
run_name = sys.argv[1]
measurements_dir="measurements"

n_cells = np.array([4, 8, 16, 32, 64, 128])

measurements_dict = {}
c_invs = []
c_sigmas = []



for filename in os.listdir(run_name):
    f = os.path.join(run_name, filename)
    # checking if it is a file
    if os.path.isfile(f):
        print(filename)
    
    #accuire info about measurement
    filename_split = filename.split("_")

    c_inv = filename_split[0]
    c_sigma = filename_split[1]
    cells = filename_split[2]
    error_kind = filename_split[3][:-4]

    error = 100.0

    #read error from file
    with open(f) as myfile:
        first_line = myfile.readline()
        print(first_line)
        error = float(first_line)

    new_tuple = (float(c_inv), float(c_sigma), int(cells), error_kind)
    measurements_dict[new_tuple] = error


error_found = measurements_dict.get((0.1, 15.0, 16, "L2"))

print(measurements_dict)
print(error_found)