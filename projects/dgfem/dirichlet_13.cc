#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>


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


void write_error_file(std::string run_name, double c_inv, int c_sigma, int num_cells, std::string error_type, double error){
    //error file
    std::setprecision(17);
    auto c_inv_str = std::to_string(c_inv);
    c_inv_str.resize(4);
    auto c_sigma_str = std::to_string(c_sigma);
    std::string out_file_name = "measurements/" + run_name + "/" + std::to_string(num_cells) + "_" + c_inv_str
                                 + "_" + c_sigma_str + "_" + error_type + ".txt";
    std::ofstream out_file(out_file_name);
    out_file << error;
    out_file.close();
}


int main(int argc, char *argv[]){

std::string run_name = argv[1];
double c_inv = std::stod(argv[2]);
double c_sigma = std::stod(argv[3]);


std::cout << "C_inv: " << c_inv << " and C_sigma: " << c_sigma << "\n"; 

//----------------------PREPARE COEFFICIENTS------------------------
// Scalar valued reaction coefficient c
auto c_coeff_lambda = [](Eigen::Vector2d x) -> double {
    return (1 + x[0]) * (1 + x[1]) * (1 + x[1]);
};
lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{c_coeff_lambda};
//Vector valued advection coefficient b
auto b_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (Eigen::Vector2d{ 2-(x[1]* x[1]) , 2 - x[0]});
};
lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{b_coeff_lambda};
// 2x2 diffusion tensor A(x)
auto a_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    double entry = std::exp(-20 * std::sqrt(x[0] * x[0] + x[1] * x[1]));
    return (Eigen::Matrix<double, 2, 2>() << entry, 0.0, 0.0, entry).finished();
};
lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{a_coeff_lambda};
//----------------------END PREPARE COEFFICIENTS------------------------

//----------------------PREPARE PRESCRIBED FUNCTIONS------------------------
// Scalar valued prescribed function gD
auto gD_lambda = [](Eigen::Vector2d x) -> double {
    return 1.0 + std::sin(M_PI * (1.0 + x[0]) * (1.0 + x[1]) * (1.0 + x[1]) * 0.125);
};
lf::dgfe::MeshFunctionGlobalDGFE m_gD{gD_lambda};

auto gN_lambda = [](Eigen::Vector2d x) -> double {
    return 0.0;
};
lf::dgfe::MeshFunctionGlobalDGFE m_gN{gN_lambda};

// Scalar valued prescribed function f
auto f_lambda = [](Eigen::Vector2d x) -> double {
    return (-0.7853981633974483*(1 + x[0])*std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)))/std::pow(E,20.0*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2))) + 
   0.7853981633974483*(2 - x[0])*(1 + x[0])*(1 + x[1])*std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)) + 
   0.39269908169872414*std::pow(1 + x[1],2)*(2 - std::pow(x[1],2))*std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)) + 
   (15.707963267948966*(1 + x[0])*x[1]*(1 + x[1])*std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)))/(std::pow(E,20.0*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2)))*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2))) + 
   (7.853981633974483*x[0]*std::pow(1 + x[1],2)*std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)))/(std::pow(E,20.0*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2)))*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2))) + 
   (0.6168502750680849*std::pow(1 + x[0],2)*std::pow(1 + x[1],2)*sdt(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)))/std::pow(E,20.0*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2))) + 
   (0.15421256876702122*std::pow(1 + x[1],4)*std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)))/std::pow(E,20.0*std::sqrt(std::pow(x[0],2) + std::pow(x[1],2))) + 
   (1 + x[0])*std::pow(1 + x[1],2)*(1 + std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)));
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
    auto l2_projection = lf::dgfe::L2ProjectionSqrtANablaBasis<double>(dgfe_space_ptr, m_a_coeff, 10);

    //run it
    auto l2_error = run_convergence(c_inv, c_sigma, 20, run_name, dgfe_space_ptr, l2_projection, m_a_coeff, m_b_coeff, m_c_coeff, m_gD, m_gN, m_f, m_gD);

    //error file
    auto n_cells = mesh_ptr->NumEntities(0);
    write_error_file(run_name, c_inv, c_sigma, n_cells, "L2", l2_error);

}



return 0;
} //end main