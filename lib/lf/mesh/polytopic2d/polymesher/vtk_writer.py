#!/usr/bin/env python
# coding: utf-8

# This code reads polytopic mesh data written by the matlab program PolyMesher: https://link.springer.com/article/10.1007/s00158-011-0706-z.
# It creates a pyvista.PolyData object out of it and writes it to vtk format.
# To prevent the meshes from changing every time the matlab script is run, append "\_final" to the matlab output files' names. 


import numpy as np
import pyvista as pv

matlab_output_folder = 'matlab_out'
output_file_name = "meshes/Mesh.vtk"

def make_mesh(matlab_out_folder, nElements):
    #read data
    nodes_raw = np.loadtxt(matlab_out_folder + '/' + 'unit_square_' + str(nElements) + '_nodes_final.txt', delimiter = ',')
    elements_str = np.loadtxt(matlab_out_folder + '/' + 'unit_square_' + str(nElements) + '_elements_final.txt', delimiter ='\n', dtype=str)
    
    #extract cell information and put into right format
    subtract_one = lambda x: x-1
    elements_list = []
    for array in elements_str:
        numpy_array = np.fromstring(array, sep=',',dtype='i')
        numpy_array = subtract_one(numpy_array)
        elements_list.append(np.insert(numpy_array, 0, numpy_array.size))
        
    cells = np.hstack(elements_list)
    #print(elements_list[57])
    #print(elements_list[68])
    #print(elements_list[438])
                             
    #add z coordinates of vertices => all 0.0
    nodes = np.column_stack((nodes_raw, np.zeros(nodes_raw.shape[0])))
    
    return pv.PolyData(nodes, cells)      

def plot_mesh(mesh):
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color = 'blue')
    plotter.camera_position = 'xy'

    dargs = dict(show_edges=True)
    plotter.add_mesh(mesh, color="#aedfc0", **dargs)
    plotter.set_background('white')

    plotter.show()


mesh = make_mesh(matlab_output_folder, 1000)
plot_mesh(mesh)
mesh.save(output_file_name, binary=False)

