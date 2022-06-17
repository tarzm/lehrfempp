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

    //setup dgfe space
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 1));
    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> sol_mesh_function(dgfe_space_ptr, sol_vec);
    
    //compute l2 error
    double l2_barycenter_error = lf::dgfe::L2ErrorBarycenter(sol_mesh_function, solution_lambda);

    std::cout << "The error of the solution via MeshFunction is: " << l2_barycenter_error << "\n";
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
    lf::dgfe::DGFEO2LocalLoadVector vectorProvider(polynomial);

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

    //setup dgfe space
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> sol_mesh_function(dgfe_space_ptr, sol_vec);
    
    //compute l2 error
    double l2_barycenter_error = lf::dgfe::L2ErrorBarycenter(sol_mesh_function, solution_lambda);

    std::cout << "The error of the solution via MeshFunction is: " << l2_barycenter_error << "\n";
}

// TEST(dgfe_01_providers, meshFunction){
//     //SIMPLE TEST MESH ----------------------------------
//     auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);


//     //setup of dofhandler
//     std::map<lf::base::RefEl, base::size_type> dofmap;
//     dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 4));
//     lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

//     //setup dgfe space
//     std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 1));


//     //galerkin matrix initialization
//     lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
//     M.setZero();

//     //initialization of element matrix provider
//     lf::dgfe::DGFEO1MassElementMatrix massMatrixProvider;

//     //assemble mass matrix
//     lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

//     //initialize polynomial for the rhs load vector: 2 + 1.5(x^2 * y)
//     lf::dgfe::Polynomial polynomial;
//     polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
//     polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
//     //initialize load vector provider
//     DGFEO1LocalLoadVector vectorProvider(polynomial);

//     //assemble rhs vector
//     Eigen::VectorXd rhs(dofhandler.NumDofs());
//     rhs.setZero();
//     lf::assemble::AssembleVectorLocally(0, dofhandler, vectorProvider, rhs);

//     //solve LSE
//     Eigen::SparseMatrix<double> M_crs = M.makeSparse();
//     Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//     solver.compute(M_crs);
//     LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
//     Eigen::VectorXd sol_vec = solver.solve(rhs);
//     LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

//     // solution lambda
//     auto solution_lambda = [](Eigen::Vector2d x) -> double {
//         return (2.0 + 1.5 * x[0] * x[0] * x[1]);
//     };

//     lf::dgfe::MeshFunctionDGFE<double> sol_mesh_function(dgfe_space_ptr, sol_vec);

//     double l2_barycenter_error = lf::dgfe::L2ErrorBarycenter(sol_mesh_function, solution_lambda);


    
//     std::cout << "The error of the solution via MeshFunction is: " << l2_barycenter_error << "\n";

// }

}