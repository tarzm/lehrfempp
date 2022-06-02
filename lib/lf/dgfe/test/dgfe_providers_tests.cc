/**
 * @file
 * @brief Tests of the the discontinuous galerkin finite elements method's assembly alogrithms
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



TEST(dgfe_O1_providers, O1massMatrixAnd01LocalLoadVector){

    //MESH FROM FILE -----------------------------
    // std::filesystem::path here = __FILE__;
    // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_100_cells.vtk";
    // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    // auto mesh_ptr = reader.mesh();

    //SIMPLE TEST MESH ----------------------------------
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //UNIT SQUARE SINGLE POLYGON MESH--------------------
    // using coord_t = Eigen::Vector2d;
    // using size_type = lf::mesh::Mesh::size_type;
    // double scale = 1.0;
    // std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
    // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
    // mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3}}, nullptr);
    // auto mesh_ptr = mesh_factory_ptr->Build();


    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 4));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();

    //initialization of element matrix provider
    lf::dgfe::DGFEO1MassElementMatrix massMatrixProvider;

    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //initialize polynomial for the rhs load vector: 2 + 1.5(x^2 * y)
    lf::dgfe::Polynomial polynomial;
    polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
    polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
    //initialize load vector provider
    DGFEO1LocalLoadVector vectorProvider(polynomial);

    //assemble rhs vector
    Eigen::VectorXd rhs(dofhandler.NumDofs());
    rhs.setZero();
    lf::assemble::AssembleVectorLocally(0, dofhandler, vectorProvider, rhs);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    // solution lambda
    auto solution_lambda = [](Eigen::Vector2d x) -> double {
        return (2.0 + 1.5 * x[0] * x[0] * x[1]);
    };

    // std::cout << std::left << std::setw(10) << "Expected" << std::right << std::setw(16)
    //         << "Actual" << std::setw(16) << "Difference" << std::setw(16) << "Area" << std::endl;
    // std::cout << "---------------------------------------------" << std::endl;
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
            int legendre_degree_x = i / 2;
            int legendre_degree_y = i % 2;
            dgfe_sol += sol_vec[dof_idx[i]] * legendre_polynomial(legendre_degree_x, mean.col(0)[0]) * legendre_polynomial(legendre_degree_y, mean.col(0)[1]);
        }

        // std::cout << std::left << std::setw(10) << true_sol << std::right << std::setw(16)
        //       << dgfe_sol << std::setw(16) << std::abs(true_sol - dgfe_sol)<< std::setw(16) << integrate(corners, 0, 0) << std::endl;

        
        //add [difference * area] to the error
        error_sum += std::abs(true_sol - dgfe_sol) * integrate(corners, 0, 0);

        //compute value at  barycenter
    }


    std::cout << "---------------------------------------------" << std::endl;

    // for (int i = 0; i < dofhandler.NumDofs(); i++){
    //     std::cout << "Coefficient " << i << " has value " << sol_vec[i] << "\n";
    // }
    std::cout << "The error of the solution is: " << error_sum << "\n";
}

TEST(dgfe_O2_providers, O2massMatrixAnd02LocalLoadVector){

    //MESH FROM FILE -----------------------------
    // std::filesystem::path here = __FILE__;
    // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_100_cells.vtk";
    // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    // auto mesh_ptr = reader.mesh();

    //SIMPLE TEST MESH ----------------------------------
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //UNIT SQUARE SINGLE POLYGON MESH--------------------
    // using coord_t = Eigen::Vector2d;
    // using size_type = lf::mesh::Mesh::size_type;
    // double scale = 1.0;
    // std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
    // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
    // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
    // mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3}}, nullptr);
    // auto mesh_ptr = mesh_factory_ptr->Build();


    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 9));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();

    //initialization of element matrix provider
    lf::dgfe::DGFEO2MassElementMatrix massMatrixProvider;

    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //initialize polynomial for the rhs load vector: 2 + 1.5(x^2 * y)
    lf::dgfe::Polynomial polynomial;
    polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
    polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
    //initialize load vector provider
    DGFEO2LocalLoadVector vectorProvider(polynomial);

    //assemble rhs vector
    Eigen::VectorXd rhs(dofhandler.NumDofs());
    rhs.setZero();
    lf::assemble::AssembleVectorLocally(0, dofhandler, vectorProvider, rhs);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    // solution lambda
    auto solution_lambda = [](Eigen::Vector2d x) -> double {
        return (2.0 + 1.5 * x[0] * x[0] * x[1]);
    };

    // std::cout << std::left << std::setw(10) << "Expected" << std::right << std::setw(16)
    //         << "Actual" << std::setw(16) << "Difference" << std::setw(16) << "Area" << std::endl;
    // std::cout << "---------------------------------------------" << std::endl;
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
        for (int i = 0; i < 9; i++){
            auto degrees = multiIndexToDegree(i , 2);
            auto legendre_degree_x = degrees.first;
            auto legendre_degree_y = degrees.second;
            for (int degree_x = 0; degree_x <= legendre_degree_x; degree_x++){
                for (int degree_y = 0; degree_y <= legendre_degree_y; degree_y++){
                    //check wether coefficients are not 0
                    scalar_t coeff_x = legendre_coeffs_(legendre_degree_x, degree_x);
                    scalar_t coeff_y = legendre_coeffs_(legendre_degree_y, degree_y);
                    if (coeff_x != 0.0 && coeff_y != 0.0){
                        dgfe_sol += sol_vec[dof_idx[i]] * coeff_x * coeff_y * legendre_polynomial(degree_x, mean.col(0)[0]) * legendre_polynomial(degree_y, mean.col(0)[1]);
                    }
                }
            }
        }

        // std::cout << std::left << std::setw(10) << true_sol << std::right << std::setw(16)
        //       << dgfe_sol << std::setw(16) << std::abs(true_sol - dgfe_sol)<< std::setw(16) << integrate(corners, 0, 0) << std::endl;

        
        //add [difference * area] to the error
        error_sum += std::abs(true_sol - dgfe_sol) * integrate(corners, 0, 0);

        //compute value at  barycenter
    }


    std::cout << "---------------------------------------------" << std::endl;

    // for (int i = 0; i < dofhandler.NumDofs(); i++){
    //     std::cout << "Coefficient " << i << " has value " << sol_vec[i] << "\n";
    // }


    std::cout << "The error of the solution is: " << error_sum << "\n";
}

}