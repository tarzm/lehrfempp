/**
 * @file
 * @brief functionalities to run convergence studies
 * @author Tarzis Maurer
 * @date August 22
 * @copyright ETH Zurich
*/

#ifndef RUN_CONVERGENCE_h
#define RUN_CONVERGENCE_h

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
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

using l2_proj_sqrt_a_nabla_basis = std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<double>>, std::vector<lf::dgfe::MeshFunctionDGFE<double>>>;


template<typename MESHFUNC_A, typename MESHFUNC_B, typename MESHFUNC_C, typename MESHFUNC_gN, typename MESHFUNC_gD,  typename MESHFUNC_f, typename MESHFUNC_true>
double run_convergence(double c_inv, double c_sigma, unsigned integration_degree, std::string run_name, std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, l2_proj_sqrt_a_nabla_basis l2_projection,
                     MESHFUNC_A &m_a, MESHFUNC_B &m_b, MESHFUNC_C &m_c, MESHFUNC_gD &m_gD, MESHFUNC_gN &m_gN, MESHFUNC_f &m_f, MESHFUNC_true &m_true){

auto n_cells = dgfe_space_ptr->Mesh()->NumEntities(0);
auto mesh_ptr = dgfe_space_ptr->Mesh();



//----------------------PREPARE BOUNDARY EDGE SETS------------------------
auto boundary_edge = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);
//boundary_N_edge
lf::mesh::utils::CodimMeshDataSet<bool> boundary_n_edge(mesh_ptr, 1, false);
//boundary_0_edge
lf::mesh::utils::CodimMeshDataSet<bool> boundary_0_edge(mesh_ptr, 1, false);
//boundary_D_edge is along x axis
lf::mesh::utils::CodimMeshDataSet<bool> boundary_d_edge(mesh_ptr, 1, false);
//boundary_minus_edge
lf::mesh::utils::CodimMeshDataSet<bool> boundary_minus_edge(mesh_ptr, 1, false);
//boundary_plus_edge
lf::mesh::utils::CodimMeshDataSet<bool> boundary_plus_edge(mesh_ptr, 1, false);

//assemble boundary_0_edge
for (auto edge : mesh_ptr->Entities(1)){
    if (boundary_edge(*edge)){
        auto corners = lf::geometry::Corners(*(edge->Geometry()));
        auto normal = lf::dgfe::outwardNormal(corners);
        auto avg_corners = corners.rowwise().mean();
        auto a = m_a.lambda()(avg_corners);
        if (normal.dot(a * normal) > 0.0){
            boundary_0_edge(*edge) = true;
        }
    }
}
//assemble minus and plus boundary
for (auto edge : mesh_ptr->Entities(1)){
    if (boundary_edge(*edge) && !boundary_0_edge(*edge)){
        //normal n
        auto corners = lf::geometry::Corners(*(edge->Geometry()));
        auto normal = lf::dgfe::outwardNormal(corners);
        auto avg_corners = corners.rowwise().mean();
        auto b = m_b.lambda()(avg_corners);
        if(b.dot(normal) < 0){
            boundary_minus_edge(*edge) = true;
        } else {
            boundary_plus_edge(*edge) = true;
        }
    }
}
//assemble boundary d and boundary n
for (auto edge : mesh_ptr->Entities(1)){
    // if (boundary_0_edge(*edge)){
    //     auto corners = lf::geometry::Corners(*(edge->Geometry()));
    //     if (corners(0,0) == 1.0 && corners(0,1) == 1.0){
    //         boundary_n_edge(*edge) = true;
    //     } else {
    //         boundary_d_edge(*edge) = true;
    //     }
    // }
    if (boundary_edge(*edge)){
        boundary_d_edge(*edge) = true;
    }
}

// std::cout << "PART OF BOUNDARY 0:\n";
// for (auto edge : mesh_ptr->Entities(1)){
//     if(boundary_0_edge(*edge)){
//         std::cout << mesh_ptr->Index(*edge) << " ";
//     }
// }
// std::cout << "\n";

// std::cout << "PART OF BOUNDARY D:\n";
// for (auto edge : mesh_ptr->Entities(1)){
//     if(boundary_d_edge(*edge)){
//         std::cout << mesh_ptr->Index(*edge) << " ";
//     }
// }
// std::cout << "\n";

// std::cout << "PART OF BOUNDARY N:\n";
// for (auto edge : mesh_ptr->Entities(1)){
//     if(boundary_n_edge(*edge)){
//         std::cout << mesh_ptr->Index(*edge) << " ";
//     }
// }
// std::cout << "\n";

// std::cout << "PART OF BOUNDARY minus:\n";
// for (auto edge : mesh_ptr->Entities(1)){
//     if(boundary_minus_edge(*edge)){
//         std::cout << mesh_ptr->Index(*edge) << " ";
//     }
// }
// std::cout << "\n";

// std::cout << "PART OF BOUNDARY plus:\n";
// for (auto edge : mesh_ptr->Entities(1)){
//     if(boundary_plus_edge(*edge)){
//         std::cout << mesh_ptr->Index(*edge) << " ";
//     }
// }
// std::cout << "\n";

    
//----------------------END PREPARE BOUNDARY EDGE SETS------------------------


//----------------------ASSEMBLE GALERKIN MATRIX------------------------
lf::dgfe::DiscontinuityPenalization disc_pen(dgfe_space_ptr, c_inv, c_sigma);
unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
//initialization of advection reaction element matrix provider
lf::dgfe::AdvectionReactionElementMatrixProvider<double, decltype(m_b), decltype(m_c), decltype(boundary_edge)>
                advectionReactionProvider(dgfe_space_ptr, m_b, m_c, boundary_edge, boundary_d_edge, boundary_minus_edge, integration_degree);

//galerkin matrix initialization
lf::assemble::COOMatrix<double> A(n_dofs, n_dofs);
A.setZero();

//assembler initialization
lf::dgfe::DiffusionMatrixAssembler<decltype(A), double, decltype(m_a), decltype(boundary_edge)>
                diffusionAssembler(dgfe_space_ptr, m_a, boundary_edge, boundary_d_edge, integration_degree, disc_pen, l2_projection);

//assemble galerkin matrix
//lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), advectionReactionProvider, A);
diffusionAssembler.assemble(A);



//----------------------END ASSEMBLE GALERKIN MATRIX------------------------


//----------------------ASSEMBLE RHS------------------------
//rhs initialization
Eigen::VectorXd rhs(n_dofs);
rhs.setZero();

//initialization of element vector provider
//lf::dgfe::AdvectionReactionDiffusionRHS<double, decltype(m_a), decltype(m_b), decltype(boundary_edge), decltype(m_f), decltype(m_gD), decltype(m_gN)>
                        //advectionReactionDiffusionRHS(dgfe_space_ptr, m_f, m_gD, m_gN, m_a, m_b, boundary_minus_edge, boundary_d_edge, boundary_n_edge, integration_degree, disc_pen, l2_projection);

//initialization of RHS Assembler
lf::dgfe::AdvectionReactionDiffusionRHSAssembler<double, decltype(m_a), decltype(m_b), decltype(boundary_edge), decltype(m_f), decltype(m_gD), decltype(m_gN), decltype(rhs)>
                        RHSAssembler(dgfe_space_ptr, m_f, m_gD, m_gN, m_a, m_b, boundary_minus_edge,
                        boundary_d_edge, boundary_n_edge, integration_degree, disc_pen, l2_projection);

//assemble rhs vector
//lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), advectionReactionDiffusionRHS, rhs);
RHSAssembler.assemble(rhs);
//----------------------END ASSEMBLE RHS------------------------


//----------------------SOLVE LSE------------------------
Eigen::SparseMatrix<double> A_crs = A.makeSparse();




Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
solver.compute(A_crs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
Eigen::VectorXd sol_vec = solver.solve(rhs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");



if (A_crs.rows() == 32 || A_crs.rows() == 8){
    //Show Galerkin Matrix
    //get dense matrix
    auto dense_A = Eigen::MatrixXd(A_crs);
    std::cout << "\n\n \t \t Galerkin Matrix: \n " << dense_A << "\n\n";

    //Show RHS vec
    //get dense matrix
    std::cout << "\n\n \t \t RHS Vec: \n " << rhs << "\n\n";


    //Show solution dof vector
    std::cout << "\n\n \t \t SOL Vec: \n " << sol_vec << "\n\n";
}
//----------------------END SOLVE LSE------------------------


//----------------------MESH FUNCTION AND ERROR CALCULATION------------------------
lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, sol_vec);

//calculate with mesh function error function
double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(dgfe_mesh_function), decltype(m_true)>(dgfe_mesh_function, m_true, mesh_ptr, 30);

// std::cout << "\n";
// std::cout << "Mesh Function error: " << mesh_func_l2_error;
// std::cout << " with C_inv: " << c_inv << " and C_sigma: " << c_sigma << " run with " << mesh_ptr->NumEntities(0) << " cells" << "\n\n";
//----------------------END MESH FUNCTION AND ERROR CALCULATION------------------------

//----------------------Show Mesh info------------------------

if (false){
    int counter = 0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        std::cout << "Cell " << counter << " has coordinates \n" << corners << "\n";
        counter++;
    }

    std::cout << "NORMALS:\n\n";

    for (auto cell : mesh_ptr->Entities(0)){
        std::cout << "CELL " << mesh_ptr->Index(*cell) << ":\n";
        int edge_count = 0;
        for (auto edge : cell->SubEntities(1)){
            //normal n
            auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
            //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
            normal *= (int) (cell->RelativeOrientations()[edge_count]);
            std::cout << "Edge " << edge_count << " has coordinates and normal \n";
            std::cout << lf::geometry::Corners(*(edge->Geometry())) << " and\n" << normal << "\n";
            edge_count++;
        }
    }


    std::cout << "MAPPINGS:\n\n";

    for (auto cell : mesh_ptr->Entities(0)){
        std::cout << "CELL " << mesh_ptr->Index(*cell) << ":\n";

        lf::dgfe::BoundingBox box(*cell);
        Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(2,1);
        reference(0,0) = .8125;
        reference(1,0) = .875;

        auto mapped = box.inverseMap(reference);
        std::cout << "CELL " << mesh_ptr->Index(*cell) << "\n";
        std::cout << "mapping maps to \n";
        std::cout << mapped << "\n";
    }
}
//----------------------END Show Mesh info------------------------


//----------------------PLOT FUNCTIONS------------------------
//output solution mesh function and true sol.
auto mesh_writer = lf::io::NumpyPolytopicWriter(mesh_ptr);
std::string meshfunc_file = "functions/problem_solution" + std::to_string(n_cells) + ".txt";
std::string truefunc_file = "functions/true_solution" + std::to_string(n_cells) + ".txt";
mesh_writer.writeSimple<decltype(dgfe_mesh_function)>(dgfe_mesh_function, meshfunc_file);
mesh_writer.writeSimple<decltype(m_true)>(m_true, truefunc_file);
//----------------------END PLOT FUNCTIONS------------------------





//----------------------WRITE ERROR TO FILE------------------------
// std::ostringstream c_inv_stream;
// c_inv_stream << std::fixed;
// c_inv_stream << std::setprecision(2);
// c_inv_stream << c_inv;
// std::string c_inv_string = c_inv_stream.str();

// std::ostringstream c_sigma_stream;
// c_sigma_stream << std::fixed;
// c_sigma_stream << std::setprecision(2);
// c_sigma_stream << c_sigma;
// std::string c_sigma_string = c_sigma_stream.str();

// std::string output_file = "measurements/" + run_name + "/" + c_inv_string + "_" + c_sigma_string + "_" + std::to_string(n_cells) + "_L2.txt";
// std::ofstream file (output_file);
// if (file.is_open()){
//     file << std::to_string(mesh_func_l2_error);
//     file.close();
// } else {
//     std::cout << "Unable to open file \n";
// }
//----------------------END WRITE ERROR TO FILE------------------------
return mesh_func_l2_error;

} //end templated function

#endif //RUN_CONVERGENCE_h