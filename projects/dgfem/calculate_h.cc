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

// Eigen::VectorXd some_vec(2);
// some_vec << 4.0, 9.0;
// Eigen::Vector2d sqrt_vec = sqrt(some_vec);
// std::cout << "Sqrt vec: " << sqrt_vec[0] << " , " << sqrt_vec[1] << "\n";

// Eigen::Matrix2d a(2,2);
// a << 1,0,0,1;

// for (int i = 0; i < 1000; i++){
//     const Eigen::Vector2d n = Eigen::Vector2d::Random();
//     auto min_n = (-1) * n;

//     double AF_p = (a * n).lpNorm<Eigen::Infinity>();
//     double AF_m = (a * min_n).lpNorm<Eigen::Infinity>();

    

//     if (!(AF_p * AF_p == AF_m * AF_m)){
//         std::cout << "One down: " << AF_p * AF_p << " and " << AF_m * AF_m << "\n";
//         std::cout << "where plus:\n " << n << "\n";
//         std::cout << "and minus:\n " << min_n << "\n\n";
//     }

// }



//prepare list of number of cells
std::vector<int> n_cells_list;
for (int i = 2; i < 13; i++){
    n_cells_list.push_back(std::pow(2,i));
}

//ouput of number cof cells in mesh
for (auto n : n_cells_list){
    std::cout << n << ",";
}
std::cout << "\n";

//loop over meshes
for (int i : n_cells_list){

    std::string num_cells = std::to_string(i);

    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_" + num_cells + "_cells.vtk";
    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    auto mesh_ptr = reader.mesh();




    //calculate h

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