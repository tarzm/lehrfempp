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
#include "run_convergence.h"

int main(int argc, char *argv[]){

std::string run_name = argv[1];
double c_inv = std::stod(argv[2]);
double c_sigma = std::stod(argv[3]);


std::cout << "C_inv: " << c_inv << " and C_sigma: " << c_sigma << "\n"; 

//----------------------PREPARE COEFFICIENTS------------------------
// Scalar valued reaction coefficient c
auto c_coeff_lambda = [](Eigen::Vector2d x) -> double {
    return 0.0;
};
lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{c_coeff_lambda};
//Vector valued advection coefficient b
auto b_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (Eigen::Vector2d{0.0, 0.0});
};
lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{b_coeff_lambda};
// 2x2 diffusion tensor A(x)
auto a_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return (Eigen::Matrix<double, 2, 2>() << 1.0, 0.0, 0.0, 1.0).finished();
};
lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{a_coeff_lambda};
//----------------------END PREPARE COEFFICIENTS------------------------

//----------------------PREPARE PRESCRIBED FUNCTIONS------------------------
// Scalar valued prescribed function gD
auto gD_lambda = [](Eigen::Vector2d x) -> double {
    return 1 + x[0] * x[0] + 2 * x[1] * x[1];
};
lf::dgfe::MeshFunctionGlobalDGFE m_gD{gD_lambda};

auto gN_lambda = [](Eigen::Vector2d x) -> double {
    return 1 + x[0] * x[0] + 2 * x[1] * x[1];
};
lf::dgfe::MeshFunctionGlobalDGFE m_gN{gN_lambda};

// Scalar valued prescribed function f
auto f_lambda = [](Eigen::Vector2d x) -> double {
    return -6.0;
};
lf::dgfe::MeshFunctionGlobalDGFE m_f{f_lambda};
//----------------------END PREPARE PRESCRIBED FUNCTIONS------------------------



//-----------------------RUN  IT------------------------------------------
//loop over meshes
for (int i = 4; i < argc; i++){

    std::string num_cells = argv[i];

    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_" + num_cells + "_cells.vtk";
    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    auto mesh_ptr = reader.mesh();

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    //Setup l2 projection of sqrt(A) * nabla(basis)
    auto l2_projection = lf::dgfe::L2ProjectionSqrtANablaBasis<double>(dgfe_space_ptr, m_a_coeff, 20);

    //run it
    run_convergence(c_inv, c_sigma, 10, run_name, dgfe_space_ptr, l2_projection, m_a_coeff, m_b_coeff, m_c_coeff, m_gD, m_gN, m_f, m_gD);

}





return 0;
} //end main