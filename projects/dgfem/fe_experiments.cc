#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <lf/mesh/mesh.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/base/base.h>
#include <lf/dgfe/dgfe.h>
#include <lf/fe/fe.h>
#include <lf/uscalfe/uscalfe.h>

#include "lf/mesh/test_utils/test_meshes.h"

int main(int /*argc*/, char ** /*argv*/) {

// abbreviations for types
using size_type = lf::base::size_type;
using glb_idx_t = lf::assemble::glb_idx_t;
using coord_t = Eigen::Vector2d;

auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3,1);

// solution function
auto alpha = [](Eigen::Vector2d x) -> double {
    return exp(x[0] * x[1]);
};
// Wrap  coefficient into a MeshFunction
lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};

// mass matrix load
auto gamma = [](Eigen::Vector2d x) -> double {
    return 1.0;
};
// Wrap  coefficient into a MeshFunction
lf::mesh::utils::MeshFunctionGlobal mf_gamma{gamma};


auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

// Dimension of finite element space`
const size_type N_dofs(dofh.NumDofs());

// Matrix in triplet format holding Galerkin matrix, zero initially.
lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
A.setZero();

//mass matrix provider
lf::fe::MassElementMatrixProvider<double, decltype(mf_gamma)> mass_matrix_provider(fe_space, mf_gamma);

//load vector provider
lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_alpha)> load_vec_provider(mf_alpha);

Eigen::VectorXd rhs(N_dofs);
rhs.setZero();

lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_matrix_provider, A);
lf::assemble::AssembleVectorLocally(0, dofh, load_vec_provider, rhs);

Eigen::SparseMatrix<double> A_crs = A.makeSparse();
// Solve linear system using Eigen's sparse direct elimination
// Examine return status of solver in case the matrix is singular
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
solver.compute(A_crs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
Eigen::VectorXd sol_vec = solver.solve(rhs);
LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);

double L2err = std::sqrt(lf::fe::IntegrateMeshFunction(*mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_alpha), 2));

std::cout << "L2 error is " << L2err << "\n";




return 0;

}