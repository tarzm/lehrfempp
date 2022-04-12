This is a short instruction to the generation of polytopic meshes and importing them to LehrFem++ with the given infrastructure.
------------------------------------------------------------------------------------------------------------------

1. ########## Mesh creation with matlab

Run *mesh.m* in matlab. The third argument of the function *write_mesh()* is the number of polytopic cells that will be present in the mesh.
This will create 2 files in the folder *matlab_out/*. One with information of the points and one with the information of the
polygons in the mesh.


2. ########## Rename files written by matlab

To prevent the created files and meshes from unwanted changes whenever mesh.m is run, please rename the files generated in *matlab_out/*.
Just add "_final" to the filenames => "unit_square_1000_elements.txt" becomes "unit_square_1000_elements_final.txt".


3. ########## Run the jupyter notebook 

Run the jupyter notebook *vtk_writer.ipynb* with the right parameters. This creates a mesh in the .vtk file format.
It should be quite self-explanatory but here the commdands to work with the existing functions to go on:


mesh = make_mesh(matlab_output_folder, 1000) #second argument is the 1000 cells again
plot_mesh(mesh)
mesh.save(output_file_name, binary=False)


4. ########## Import the mesh into LehrFEM++

Read the mesh from the .vtk file with LehrFEM++ functionalities:


std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
std::filesystem::path here = __FILE__;
auto mesh_file = (here.parent_path() / "YOUR_RELATIVE_PATH_TO_THE_MESH").string();
VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
auto mesh = reader.mesh();

	




