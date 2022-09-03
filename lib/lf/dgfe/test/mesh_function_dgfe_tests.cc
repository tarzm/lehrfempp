/**
 * @file
 * @brief Tests of the dgfe mesh function functionalities
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
 */

#include <cmath>
#include <filesystem>
#include <limits>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include <lf/assemble/assemble.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(meshFunction, singleSquareO2L2ErrorSubTessellation){

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

    //setup dof vector
    Eigen::VectorXd dof_vector = Eigen::VectorXd::Zero(9);
    dof_vector[0] = 0.5;
    dof_vector[3] = 0.5;

    //true solution is x
    auto true_sol_lambda = [](Eigen::Vector2d x) -> double {
        return x[0];
    };

    //setup dgfe space with legendre polynomials of degree 2
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));

    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vector);
    //setup global mesh function
    lf::dgfe::MeshFunctionGlobalDGFE<decltype(true_sol_lambda)> lambda_mesh_func(true_sol_lambda);

    //calculate with mesh function error function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double>(dgfe_mesh_function, lambda_mesh_func, mesh_ptr, 2);

    //mesh function should be exact
    EXPECT_NEAR(0.0, mesh_func_l2_error, std::numeric_limits<double>::epsilon());
}

TEST(meshFunctionGrad, singleSquareO2L2ErrorSubTessellation){

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

    //setup dof vector
    Eigen::VectorXd dof_vector = Eigen::VectorXd::Zero(9);
    dof_vector[0] = 0.5;
    dof_vector[3] = 0.5;

    //true solution is x
    auto true_sol_lambda = [](Eigen::Vector2d x) -> double {
        return x[0];
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_true{true_sol_lambda};

    //true grad function is [1,0]^T
    auto true_sol_grad_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d() << 1.0, 0.0).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_true_grad{true_sol_grad_lambda};

    //setup dgfe space with legendre polynomials of degree 2
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));

    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vector);
    //setup dgfe mesh gradient function
    lf::dgfe::MeshFunctionGradDGFE<double> dgfe_mesh_function_grad(dgfe_space_ptr, dof_vector);
    //setup global mesh function
    lf::dgfe::MeshFunctionGlobalDGFE<decltype(true_sol_lambda)> lambda_mesh_func(true_sol_lambda);

    //calculate with mesh function error function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double>(dgfe_mesh_function, lambda_mesh_func, mesh_ptr, 2);

    //mesh function should be exact
    EXPECT_NEAR(0.0, mesh_func_l2_error, std::numeric_limits<double>::epsilon());
}


TEST(meshFunctionGrad, LSETest){
    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_100_cells.vtk";
    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    auto mesh_ptr = reader.mesh();

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    // Exact solution u
    auto u_lambda = [](Eigen::Vector2d x) -> double {
        return std::log(x[0] * x[0] + x[1] + 1.0);
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_u{u_lambda};

    // Gradient of exact solution
    auto grad_u_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
        double den = x[0] * x[0] + x[1] + 1.0;
        return ((Eigen::Vector2d() << 2.0 * x[0], 1.0).finished()) / den;
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_grad_u{grad_u_lambda};

    auto n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();

    //mass matrix initialization
    lf::assemble::COOMatrix<double> M(n_dofs, n_dofs);
    M.setZero();

    //initialization of element matrix provider
    lf::dgfe::DGFEMassElementMatrixST<double> massMatrixProvider(10, 2);

    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), massMatrixProvider, M);

    //rhs vector initialization
    Eigen::VectorXd rhs(n_dofs);
    rhs.setZero();
    //initialization of vector provider
    lf::dgfe::DGFELoadElementVectorProvider<double, decltype(m_u)> vec_provider(dgfe_space_ptr, m_u);
    //assemble load vector
    lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), vec_provider, rhs);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    //setup mesh function with solution vector
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, sol_vec);

    //setup gradient mesh function with solution vector
    lf::dgfe::MeshFunctionGradDGFE<double> dgfe_grad_mesh_function(dgfe_space_ptr, sol_vec);

    //calculate with mesh function error function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double>(dgfe_mesh_function, m_u, mesh_ptr, 30);

    //calculate with mesh function error function
    double mesh_func_grad_l2_error = lf::dgfe::L2ErrorGradSubTessellation<double, decltype(dgfe_grad_mesh_function), decltype(m_grad_u)>(dgfe_grad_mesh_function, m_grad_u, mesh_ptr, 30);
    
    std::cout << "The error of the solution via MeshFunction is: " << mesh_func_l2_error << "\n";

    std::cout << "The error of the gradient of the solution via MeshFunction is: " << mesh_func_grad_l2_error << "\n";

}

} // namespace lf::dgfe::test