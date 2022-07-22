/**
 * @file
 * @brief Tests of the AdvectionReactionElementMatrixProvider class
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <limits>

#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include "lf/mesh/test_utils/test_meshes.h"

#include<Eigen/SparseCholesky>

namespace lf::dgfe::test {


TEST(diffusion_provider, basicEvaluation){
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    //----------------------PREPARE COEFFICIENTS------------------------
    // Scalar valued reaction coefficient c
    auto c_coeff_lambda = [](Eigen::Vector2d x) -> double {
        return 0.0;
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{c_coeff_lambda};
    //Vector valued advection coefficient b
    auto b_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d{1.0, 1.0});
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{b_coeff_lambda};
    // 2x2 diffusion tensor A(x)
    auto a_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
        return (Eigen::Matrix<double, 2, 2>() << 0.0, 0.0, 0.0, 0.0).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{a_coeff_lambda};
    //----------------------END PREPARE COEFFICIENTS------------------------

    //----------------------PREPARE BOUNDARY EDGE SETS------------------------
    auto boundary_edge = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);
    //boundary_N_edge is empty
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_n_edge(mesh_ptr, 1, false);
    //boundary_0_edge is empty
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_0_edge(mesh_ptr, 1, false);
    //boundary_D_edge is along x axis
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_d_edge(mesh_ptr, 1, false);
    //boundary_minus_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_minus_edge(mesh_ptr, 1, false);
    //boundary_plus_edge
    lf::mesh::utils::CodimMeshDataSet<bool> boundary_plus_edge(mesh_ptr, 1, false);
    for (auto edge : mesh_ptr->Entities(1)){
        if (boundary_edge(*edge)){
            //plus or minus
            //normal n
            auto corners = lf::geometry::Corners(*(edge->Geometry()));
            auto normal = lf::dgfe::outwardNormal(corners);
            auto b = b_coeff_lambda(corners.col(0));
            if(b.dot(normal) < 0){
                boundary_minus_edge(*edge) = true;
            } else {
                boundary_plus_edge(*edge) = true;
            }

            //part of boundary_d ?
            if (corners(1,0) == 0.0 && corners(1,1) == 0.0){
                boundary_d_edge(*edge) = true;
            }
        }
    }
    //----------------------END PREPARE BOUNDARY EDGE SETS------------------------

    //----------------------PREPARE PRESCRIBED FUNCTIONS------------------------
    // Scalar valued prescribed function gD
    auto gD_lambda = [](Eigen::Vector2d x) -> double {
        return std::sin(x[0] * x[0] + x[1]);
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_gD{gD_lambda};

    // Scalar valued prescribed function f
    auto f_lambda = [](Eigen::Vector2d x) -> double {
        return 1.0;
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_f{f_lambda};
    //----------------------END PREPARE PRESCRIBED FUNCTIONS------------------------

    //----------------------ASSEMBLE GALERKIN MATRIX------------------------
    double c_inv = 0.5;
    double c_sigma = 2.0;
    lf::dgfe::DiscontinuityPenalization disc_pen(dgfe_space_ptr, c_inv, c_sigma);
    unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
    //initialization of element matrix provider
    lf::dgfe::DiffusionElementMatrixProvider<double, decltype(m_a_coeff), decltype(boundary_edge)>
                    diffusionProvider(dgfe_space_ptr, m_a_coeff, boundary_edge, boundary_d_edge, 10, disc_pen);
    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> A(n_dofs, n_dofs);
    A.setZero();
    //assemble galerkin matrix
    lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), diffusionProvider, A);
    //----------------------END ASSEMBLE GALERKIN MATRIX------------------------
}

TEST(diffusion_assembler, fullDiffusionLSE){
    //UNIT SQUARE SINGLE POLYGON MESH--------------------
// using coord_t = Eigen::Vector2d;
// using size_type = lf::mesh::Mesh::size_type;
// double scale = 1.0;
// std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
// mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
// mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
// mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
// mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
// mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,4>{{0,1,2,3}}, nullptr);
// auto mesh_ptr = mesh_factory_ptr->Build();

auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

////get mesh
// std::filesystem::path here = __FILE__;
// auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_100_cells.vtk";
// lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
// auto mesh_ptr = reader.mesh();

//dgfe space
lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

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
        auto a = a_coeff_lambda(avg_corners);
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
        auto b = b_coeff_lambda(avg_corners);
        if(b.dot(normal) < 0){
            boundary_minus_edge(*edge) = true;
        } else {
            boundary_plus_edge(*edge) = true;
        }
    }
}
//assemble boundary d and boundary n
for (auto edge : mesh_ptr->Entities(1)){
    if (boundary_0_edge(*edge)){
        //normal n
        // auto corners = lf::geometry::Corners(*(edge->Geometry()));
        // auto normal = lf::dgfe::outwardNormal(corners);
        // auto avg_corners = corners.rowwise().mean();
        // auto b = b_coeff_lambda(avg_corners);
        // if(b.dot(normal) < 0){
        //     boundary_minus_edge(*edge) = true;
        // } else {
        //     boundary_plus_edge(*edge) = true;
        // }
        boundary_d_edge(*edge) = true;
    }
}

std::cout << "PART OF BOUNDARY 0:\n";
for (auto edge : mesh_ptr->Entities(1)){
    if(boundary_0_edge(*edge)){
        std::cout << mesh_ptr->Index(*edge) << " ";
    }
}
std::cout << "\n";

std::cout << "PART OF BOUNDARY D:\n";
for (auto edge : mesh_ptr->Entities(1)){
    if(boundary_d_edge(*edge)){
        std::cout << mesh_ptr->Index(*edge) << " ";
    }
}
std::cout << "\n";

std::cout << "PART OF BOUNDARY N:\n";
for (auto edge : mesh_ptr->Entities(1)){
    if(boundary_n_edge(*edge)){
        std::cout << mesh_ptr->Index(*edge) << " ";
    }
}
std::cout << "\n";

std::cout << "PART OF BOUNDARY minus:\n";
for (auto edge : mesh_ptr->Entities(1)){
    if(boundary_minus_edge(*edge)){
        std::cout << mesh_ptr->Index(*edge) << " ";
    }
}
std::cout << "\n";

std::cout << "PART OF BOUNDARY plus:\n";
for (auto edge : mesh_ptr->Entities(1)){
    if(boundary_plus_edge(*edge)){
        std::cout << mesh_ptr->Index(*edge) << " ";
    }
}
std::cout << "\n";

    
//----------------------END PREPARE BOUNDARY EDGE SETS------------------------

//----------------------PREPARE PRESCRIBED FUNCTIONS------------------------
// Scalar valued prescribed function gD
auto gD_lambda = [](Eigen::Vector2d x) -> double {
    return 1 + x[0] * x[0] + 2 * x[1] * x[1];
};
lf::dgfe::MeshFunctionGlobalDGFE m_gD{gD_lambda};

// Scalar valued prescribed function f
auto f_lambda = [](Eigen::Vector2d x) -> double {
    return -6.0;
};
lf::dgfe::MeshFunctionGlobalDGFE m_f{f_lambda};
//----------------------END PREPARE PRESCRIBED FUNCTIONS------------------------

//----------------------ASSEMBLE GALERKIN MATRIX------------------------

double c_inv = 0.15;
double c_sigma = 6.0;

double smallest_error = 20.0;
double smallest_c_inv = 10.0;
double smallest_c_sigma = 10.0;

lf::dgfe::DiscontinuityPenalization disc_pen(dgfe_space_ptr, c_inv, c_sigma);
unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
//initialization of advection reaction element matrix provider
lf::dgfe::AdvectionReactionElementMatrixProvider<double, decltype(m_b_coeff), decltype(m_c_coeff), decltype(boundary_edge)>
                advectionReactionProvider(dgfe_space_ptr, m_b_coeff, m_c_coeff, boundary_edge, boundary_d_edge, boundary_minus_edge, 10);

//galerkin matrix initialization
lf::assemble::COOMatrix<double> A(n_dofs, n_dofs);
A.setZero();
//assemble galerkin matrix

//assembler initialization
lf::dgfe::DiffusionMatrixAssembler<decltype(A), double, decltype(m_a_coeff), decltype(boundary_edge)>
                diffusionAssembler(dgfe_space_ptr, m_a_coeff, boundary_edge, boundary_d_edge, 10, disc_pen);

lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), advectionReactionProvider, A);
diffusionAssembler.assemble(A); 
//----------------------END ASSEMBLE GALERKIN MATRIX------------------------


//----------------------ASSEMBLE RHS------------------------
//initialization of element vector provider
lf::dgfe::AdvectionReactionDiffusionRHS<double, decltype(m_a_coeff), decltype(m_b_coeff), decltype(boundary_edge), decltype(m_f), decltype(m_gD), decltype(m_f)>
                        advectionReactionDiffusionRHS(dgfe_space_ptr, m_f, m_gD, m_f, m_a_coeff, m_b_coeff, boundary_minus_edge, boundary_d_edge, boundary_n_edge, 10, disc_pen);
//rhs initialization
Eigen::VectorXd rhs(n_dofs);
rhs.setZero();
//assemble rhs vector
lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), advectionReactionDiffusionRHS, rhs);
//----------------------END ASSEMBLE RHS------------------------


// //----------------------SOLVE LSE------------------------
// Eigen::SparseMatrix<double> A_crs = A.makeSparse();
// Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
// solver.compute(A_crs);
// LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
// Eigen::VectorXd sol_vec = solver.solve(rhs);
// LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
// //----------------------END SOLVE LSE------------------------

Eigen::SparseMatrix<double> A_crs = A.makeSparse();
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
solver.compute(A_crs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
Eigen::VectorXd sol_vec = solver.solve(rhs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

//----------------------MESH FUNCTION AND ERROR CALCULATION------------------------
lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, sol_vec);

//calculate with mesh function error function
double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(m_gD)>(dgfe_mesh_function, m_gD, 17);

std::cout << "Mesh Function error: " << mesh_func_l2_error;
std::cout << " with C_inv: " << c_inv << " and C_sigma: " << c_sigma << "\n\n";

if (mesh_func_l2_error < smallest_error){
    smallest_error = mesh_func_l2_error;
    smallest_c_inv = c_inv;
    smallest_c_sigma = c_sigma;
}


std::cout << "SMALLEST Mesh Function error: " << smallest_error;
        std::cout << " with C_inv: " << smallest_c_inv << " and C_sigma: " << smallest_c_sigma << "\n\n";
}


} //namespace lf::dgfe::test