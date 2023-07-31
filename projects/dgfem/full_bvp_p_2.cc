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

//function to write error to a file
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

//############ USAGE of this script in bash #########################
/*

./projects.dgfem.full_bvp <RUN_NAME> <C_INV> <C_SIGMA> <LIST OF # cells in mesh separated by space>

EXAMPLE:

./projects.dgfem.full_bvp my_run_name 0.5 20 4 8 16 32 64 128 256 512 1024 2048 4096

*/


//get run arguments: name and constants for discontinuity penalty
std::string run_name = argv[1];
double c_inv = std::stod(argv[2]);
double c_sigma = std::stod(argv[3]);
std::cout << "C_inv: " << c_inv << " and C_sigma: " << c_sigma << "\n";

//set the integration degree used for assembling Galerkin Matrix and RHS
int integration_degree = 15;

//----------------------PREPARE COEFFICIENTS------------------------
// Scalar valued reaction coefficient c
auto c_coeff_lambda = [](Eigen::Vector2d x) -> double {
    return (1 + x[0]) * (1 + x[1]) * (1 + x[1]);
};
lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{c_coeff_lambda};
//Vector valued advection coefficient b
auto b_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (Eigen::Vector2d{2 - x[0]* x[1] , 2 - x[0]*x[0]});
};
lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{b_coeff_lambda};
// Scalar valued div of advection coeff
auto div_b_coeff_lambda = [](Eigen::Vector2d x) -> double {
    return -x[1];
};
lf::dgfe::MeshFunctionGlobalDGFE m_div_b_coeff{div_b_coeff_lambda};
// 2x2 diffusion tensor A(x)
auto a_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    double entry = 1.0 + (std::sin(4.0 * M_PI *(x[0] + x[1]))) * (std::sin(4.0 * M_PI *(x[0] + x[1])));
    return (Eigen::Matrix<double, 2, 2>() << entry, 0.0, 0.0, entry).finished();
};
lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{a_coeff_lambda};
//----------------------END PREPARE COEFFICIENTS------------------------


//----------------------PREPARE PRESCRIBED FUNCTIONS------------------------
//define eulers number
const double E = 2.7182818284590452353602874713526624977572470936999595749669676277;
// Scalar valued prescribed function gD
auto gD_lambda = [](Eigen::Vector2d x) -> double {
    return 1.0 + std::sin(M_PI * (1.0 + x[0]) * (1.0 + x[1]) * (1.0 + x[1]) * 0.125);
};
lf::dgfe::MeshFunctionGlobalDGFE m_gD{gD_lambda};

auto gN_lambda = [E](Eigen::Vector2d x) -> double {
    return 0.39269908169872414*std::pow(1 + x[1],2)*
        std::cos(0.39269908169872414*(1.0 + x[0])*std::pow(1 + x[1],2))*
        (1.0 + std::pow(std::sin(4*M_PI*(x[0] + x[1])),2));
};
lf::dgfe::MeshFunctionGlobalDGFE m_gN{gN_lambda};

// Scalar valued prescribed function f
auto f_lambda = [](Eigen::Vector2d x) -> double {
    return 0.7853981633974483*(1 + x[0])*(2 - std::pow(x[0],2))*(1 + x[1])*
    std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)) + 
   0.39269908169872414*std::pow(1 + x[1],2)*(2 - x[0]*x[1])*
    std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2)) - 
   x[1]*(1 + std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))) + 
   (1 + x[0])*std::pow(1 + x[1],2)*
    (1 + std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))) - 
   19.739208802178716*(1 + x[0])*(1 + x[1])*
    std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))*
    std::cos(4*M_PI*(x[0] + x[1]))*std::sin(4*M_PI*(x[0] + x[1])) - 
   9.869604401089358*std::pow(1 + x[1],2)*
    std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))*
    std::cos(4*M_PI*(x[0] + x[1]))*std::sin(4*M_PI*(x[0] + x[1])) - 
   0.7853981633974483*(1 + x[0])*
    std::cos(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))*
    (1 + std::pow(std::sin(4*M_PI*(x[0] + x[1])),2)) + 
   0.6168502750680849*std::pow(1 + x[0],2)*std::pow(1 + x[1],2)*
    std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))*
    (1 + std::pow(std::sin(4*M_PI*(x[0] + x[1])),2)) + 
   0.15421256876702122*std::pow(1 + x[1],4)*
    std::sin(0.39269908169872414*(1 + x[0])*std::pow(1 + x[1],2))*
    (1 + std::pow(std::sin(4*M_PI*(x[0] + x[1])),2));
};
lf::dgfe::MeshFunctionGlobalDGFE m_f{f_lambda};
//----------------------END PREPARE PRESCRIBED FUNCTIONS------------------------

//-----------------------LOOP over Meshes------------------------------------------
//loop over meshes
for (int i = 4; i < argc; i++){

    //read number of cells
    std::string num_cells = argv[i];

    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_" + num_cells + "_cells.vtk";
    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    auto mesh_ptr = reader.mesh();

    //write mesh for python drawing
    //lf::io::writeMatplotlib(*mesh_ptr, "./csvs/" + std::to_string(mesh_ptr->NumEntities(0)) + ".csv");

    //dgfe space p = 2
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 1);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    //----------------------PREPARE BOUNDARY EDGE SETS------------------------
    auto boundary_edge = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);
    //boundary_N_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_n_edge(mesh_ptr, 1, false);
    //boundary_0_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_0_edge(mesh_ptr, 1, false);
    //boundary_D_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_d_edge(mesh_ptr, 1, false);
    //boundary_minus_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_minus_edge(mesh_ptr, 1, false);
    //boundary_plus_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_plus_edge(mesh_ptr, 1, false);

    //setup qr rule for segments
    const lf::quad::QuadRule qr_s = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), integration_degree);
    // qr points
    const Eigen::MatrixXd zeta_ref_s{qr_s.Points()};
    //weights
    Eigen::VectorXd w_ref_s{qr_s.Weights()};

    //BOUNDARY SETS ASSEMBLY
    for (auto cell : mesh_ptr->Entities(0)){
        for (auto edge : cell->SubEntities(1)){
            if (boundary_edge(*edge)){
                //normal n
                auto polygon_pair = dgfe_space_ptr->AdjacentPolygons(edge);
                auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
                //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
                normal *= (int) (cell->RelativeOrientations()[polygon_pair.first.second]);

                lf::dgfe::BoundingBox box(*cell);
                // qr points mapped to segment
                Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
                //gramian determinants
                Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};
                auto a_evaluated = m_a_coeff(*cell, zeta_box_s);
                double boundary_0_sum = 0.0;

                for (int i = 0; i < gram_dets_s.size(); i++){
                    boundary_0_sum += normal.dot(a_evaluated[i] * normal) * gram_dets_s[i] * w_ref_s[i];
                }
                //BOUNDARY 0 #############
                if (boundary_0_sum > 0){
                    boundary_0_edge(*edge) = true;
                    //HERE DIRICHLET AND NEUMANN ################
                    auto corners = lf::geometry::Corners(*(edge->Geometry()));
                    if(corners(0,0) == 1.0 && corners(0,1) == 1.0){ //whole edge on side x = 1
                        boundary_n_edge(*edge) = true;
                    } else {
                        boundary_d_edge(*edge) = true;
                    }
                } else { //BOUNDARY_plus and BOUNDARY_minus ###############
                    auto b_evaluated = m_b_coeff(*cell, zeta_box_s);
                    double boundary_plus_sum = 0.0;
                    for (int i = 0; i < gram_dets_s.size(); i++){
                        boundary_0_sum += b_evaluated[i].dot(normal) * gram_dets_s[i] * w_ref_s[i];
                    }
                    if (boundary_plus_sum < 0){
                        boundary_minus_edge(*edge) = true;
                    } else {
                        boundary_plus_edge(*edge) = true;
                    }
                }
            }   
        }
    }
    //----------------------END PREPARE BOUNDARY EDGE SETS------------------------


    //----------------------ASSEMBLE GALERKIN MATRIX & RHS------------------------
    //set up discontinuity penalization
    lf::dgfe::DiscontinuityPenalization disc_pen(dgfe_space_ptr, c_inv, c_sigma);
    unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> A(n_dofs, n_dofs);
    A.setZero();
    //diffusion assembler
    lf::dgfe::DiffusionMatrixAssembler<decltype(A), double, decltype(m_a_coeff), decltype(boundary_edge)>
                    diffusionAssembler(dgfe_space_ptr, m_a_coeff, boundary_edge, boundary_d_edge, integration_degree, disc_pen);
    //advection reaction matrix assembler
    lf::dgfe::AdvectionReactionMatrixAssembler<decltype(A), double, decltype(m_b_coeff), decltype(m_c_coeff), decltype(boundary_edge), decltype(m_div_b_coeff)>
                    advectionReactionAssembler(dgfe_space_ptr, m_b_coeff, m_div_b_coeff, m_c_coeff, boundary_edge, boundary_d_edge, boundary_minus_edge, integration_degree);
    //assemble matrix
    diffusionAssembler.assemble(A);
    advectionReactionAssembler.assemble(A);

    //rhs initialization
    Eigen::VectorXd rhs(n_dofs);
    rhs.setZero();
    //RHS Assembler
    lf::dgfe::AdvectionReactionDiffusionRHSAssembler<double, decltype(m_a_coeff), decltype(m_b_coeff), decltype(boundary_edge), decltype(m_f), decltype(m_gD),
                                                    decltype(m_gN), decltype(rhs)>
                            rhsAssembler(dgfe_space_ptr, m_f, m_gD, m_gN, m_a_coeff, m_b_coeff, boundary_minus_edge,
                            boundary_d_edge, boundary_n_edge, integration_degree, disc_pen);
    //assemble rhs vector
    rhsAssembler.assemble(rhs);
    //----------------------END ASSEMBLE GALERKIN MATRIX & RHS------------------------

    //----------------------SOLVE LSE------------------------
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
    //----------------------END SOLVE LSE------------------------

    //----------------------MESH FUNCTION AND ERROR CALCULATION------------------------
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, sol_vec);

    //calculate L2 error of solution with "overkill" QR
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(dgfe_mesh_function), decltype(m_gD)>(dgfe_mesh_function, m_gD, mesh_ptr, 30);
    //----------------------END MESH FUNCTION AND ERROR CALCULATION------------------------

    //----------------------WRITE ERROR TO FILE------------------------
    write_error_file(run_name, c_inv, c_sigma, std::stoi(num_cells), "L2", mesh_func_l2_error);

    std::cout << "L2 Error for " << num_cells << " cells: \t" << mesh_func_l2_error << "\n";
    //---------------------END WRITE ERROR TO FILE------------------------

}
//-----------------------END LOOP over Meshes------------------------------------------


return 0;
} //end main