#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/base/base.h>
#include <lf/dgfe/dgfe.h>
#include <lf/fe/fe.h>
#include <lf/uscalfe/uscalfe.h>

#include "lf/mesh/test_utils/test_meshes.h"
#include <Eigen/Eigen>


int main(int argc, char *argv[]){

//prepare list of number of cells
std::vector<int> n_cells_list;
for (int i = 2; i < 13; i++){
    n_cells_list.push_back(std::pow(2,i));
}

//loop over meshes
for (int i : n_cells_list){

    std::string num_cells = std::to_string(i);

    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_" + num_cells + "_cells.vtk";
    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    auto mesh_ptr = reader.mesh();

    //calculate h => biggest euclidean distance between two nodes of a polygon
    double h = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        auto n_nodes = corners.cols();
        for (int idx_0 = 0; idx_0 < n_nodes; idx_0++){
            for (int idx_1 = 0; idx_1 < n_nodes; idx_1++){

                double norm = (corners.col(idx_0)- corners.col(idx_1)).norm();
                if (norm > h){
                    h = norm;
                }
            }
        }
    }
    std::cout << h  << ", ";
}
std::cout << "\n";
return 0;
} //end main