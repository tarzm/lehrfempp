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



} //namespace lf::dgfe::test