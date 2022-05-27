/**
 * @file
 * @brief Tests of the the discontinuous galerkin finite elements
 * methods and their assembly alogrithms
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
 */

#include <cmath>
#include <filesystem>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include <lf/assemble/assemble.h>
#include "lf/mesh/test_utils/test_meshes.h"



namespace lf::dgfe::test {

using scalar_t = lf::dgfe::scalar_t;

TEST(legendre_dgfe, legendrePolynomials){
    scalar_t x = 0.8;
    scalar_t y = -0.7;
    EXPECT_EQ(legendre_polynomial(1, x) * legendre_polynomial(2,y), 0.18799999999999992);
}

TEST(legendre_dgfe, O1BasisFunctions){
    //test implementation of getting the correct exponents in the degree O1 basis functions
    //from the multi index of the basis function
    std::vector<std::pair<int, int>> check;
    check.push_back(std::make_pair(0,0));
    check.push_back(std::make_pair(0,1));
    check.push_back(std::make_pair(1,0));
    check.push_back(std::make_pair(1,1));

    for (int i = 0; i< 4; i++){
        int basis_degree_x = i / 2;
        int basis_degree_y = i % 2;
        EXPECT_EQ(std::make_pair(basis_degree_x, basis_degree_y), check[i]);
    }
}


TEST(legendre_dgfe, massMatrixAnd01LocalLoadVector){
    // //get mesh
    // std::filesystem::path here = __FILE__;
    // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_100_cells.vtk";
    // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    // auto mesh_ptr = reader.mesh();

    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);


    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 4));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();

    //initialization of element matrix provider
    lf::dgfe::DGFEO1MassElementMatrix massMatrixProvider;

    //lf::assemble::AssembleMatrixLogger()->set_level(spdlog::level::debug);
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //initialize polynomial for the rhs load vector: 2 + 1.5( x^2 * y)
    lf::dgfe::DGFEO1LocalLoadVector::Polynomial polynomial;
    polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
    polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
    //initialize load vector provider
    DGFEO1LocalLoadVector vectorProvider(polynomial);

    //assemble rhs vector
    Eigen::VectorXd rhs(dofhandler.NumDofs()) ;
    lf::assemble::AssembleVectorLocally(0, dofhandler, vectorProvider, rhs);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    std::cout << "M---------------------\n";
    std::cout << Eigen::MatrixXd(M_crs) << "\n";
    std::cout << "--------------------\n";

    std::cout << "RHS-----------------\n"; 
    for (int i = 0; i  < rhs.size(); i++){
        std::cout << rhs[i] << "\n";
    }
    std::cout << "--------------------\n";

    // solution lambda
    auto solution_lambda = [](Eigen::Vector2d x) -> double {
        return (2.0 + 1.5 * x[0] * x[0] * x[1]);
    };

    std::cout << std::left << std::setw(10) << "Expected" << std::right << std::setw(16)
            << "Actual" << std::setw(16) << "Difference" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    //compute error: I: only barycenter of cells
    double error_sum = 0;
    for (auto cell : mesh_ptr->Entities(0)){
        BoundingBox box(*cell);
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        //barycenter
        Eigen::MatrixXd mean = corners.rowwise().mean();
        
        //true solution of barycenter
        double true_sol = solution_lambda(mean.col(0));

        //calculated solution of dgfe space
        //dof indices of cell
        nonstd::span<const gdof_idx_t> dof_idx(dofhandler.GlobalDofIndices(*cell));
        double dgfe_sol = 0;
        //loop over basis functions
        for (int i = 0; i < 4; i++){
            int basis_degree_x = i / 2;
            int basis_degree_y = i % 2;
            dgfe_sol += sol_vec[dof_idx[i]] * legendre_polynomial(basis_degree_x, mean.col(0)[0]) * legendre_polynomial(basis_degree_y, mean.col(0)[1]);
        }

        std::cout << std::left << std::setw(10) << true_sol << std::right << std::setw(16)
              << dgfe_sol << std::setw(16) << std::abs(true_sol - dgfe_sol) << std::endl;

        
        //add [difference * area] to the error
        error_sum += std::abs(true_sol - dgfe_sol) * integrate(corners, 0, 0);

        //compute value at  barycenter
    }


    std::cout << "---------------------------------------------" << std::endl;

    for (int i = 0; i < dofhandler.NumDofs(); i++){
        std::cout << "Coefficient " << i << " has value " << sol_vec[i] << "\n";
    }


    std::cout << "The error of the solution is: " << error_sum << "\n";


}


} //namespace lf::dgfe:test