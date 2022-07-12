/**
 * @file
 * @brief Tests of the AdvectionReactionElementMatrixProvider class
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>

#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(advectionReaction, basicEvaluation){

//UNIT SQUARE SINGLE POLYGON MESH--------------------
using coord_t = Eigen::Vector2d;
using size_type = lf::mesh::Mesh::size_type;
double scale = 1.0;
std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,4>{{0,1,2,3}}, nullptr);
auto mesh_ptr = mesh_factory_ptr->Build();

//auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

//dgfe space
lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

// Scalar valued reaction coefficient c
auto c_coeff_lambda = [](Eigen::Vector2d x) -> double {
    return 1.0;
};
lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{c_coeff_lambda};

//Vector valued advection coefficient b
auto b_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (Eigen::Vector2d{1.0, 1.0});
};
lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{b_coeff_lambda};


auto edge_pred_0 = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);
auto edge_pred_1 = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);
auto edge_pred_2 = lf::mesh::utils::flagEntitiesOnBoundary(mesh_ptr, 1);

//initialization of element matrix provider
lf::dgfe::AdvectionReactionElementMatrixProvider<double, decltype(m_b_coeff), decltype(m_c_coeff), decltype(edge_pred_0)> advectionReactionProvider(dgfe_space_ptr, m_b_coeff, m_c_coeff, edge_pred_0, edge_pred_1, edge_pred_2, 10);

unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
//galerkin matrix initialization
lf::assemble::COOMatrix<double> A(n_dofs, n_dofs);
A.setZero();

//assemble galerkin matrix
lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), advectionReactionProvider, A);

}

} //namespace lf::dgfe::test