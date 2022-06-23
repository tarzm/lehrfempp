/**
 * @file
 * @brief Tests of the the discontinuous galerkin finite elements method's assembly alogrithms
 * @author Tarzis Maurer
 * @date May 22
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

TEST(dgfe_SubTessellation_providers, massMatrixProviderO1){
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

    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 4));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();
    //initialization of element matrix provider
    lf::dgfe::DGFEMassElementMatrixST massMatrixProvider(5, 1); //quadrule of degree 4 would be ok too but then coefficients inside the matrix would be bigger than epsilon
    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //entries of the check matrix => functions integrated over (0,1)^2
    double x = 0.5;
    double y = 0.5;
    double xy = 0.25;
    double x_2 = 1.0 / 3.0; //x^2
    double y_2 = 1.0 / 3.0; //y^2
    double x_2_y = 1.0 / 6.0; //x^2 * y
    double x_y_2 = 1.0 / 6.0; //x* y^2
    double x_2_y_2 = 1.0 / 9.0; // x^2 * y^2

    Eigen::Matrix4d check_M;
    check_M <<  1,      y,      x,      xy,
                y,      y_2,    xy,     x_y_2,
                x,      xy,     x_2,    x_2_y,
                xy,     x_y_2,  x_2_y,  x_2_y_2;

    auto M_dense = M.makeDense();


    EXPECT_TRUE(M_dense.isApprox(check_M, std::numeric_limits<double>::epsilon()));
    
}

TEST(dgfe_SubTessellation_providers, massMatrixProviderO2){
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

    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 9));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();
    //initialization of element matrix provider
    lf::dgfe::DGFEMassElementMatrixST massMatrixProvider(9, 1); 
    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    /**entries of the check matrix => functions integrated over (0,1)^2
     *BASIS FUNCTION IDX| FUNCTION EXPRESSION
     * 0                    1
     * 1                    y
     * 2                    1.5y^2 - 0.5
     * 3                    x
     * 4                    xy
     * 5                    x * (1.5y^2 - 0.5)
     * 6                    1.5x^2 - 0.5
     * 7                    y * (1.5x^2 - 0.5)
     * 8                    (1.5x^2 - 0.5) * (1.5y^2 - 0.5)
     * 
     * Now the entry (i, j) of the elements matrix is basis_i(x,y) * basis_j(x,y) integrated over (0,1)^2
    */

    double e_0_0 = 1.0;
    double e_1_0 = 0.5;
    double e_2_0 = 0.0;
    double e_3_0 = 0.5;
    double e_4_0 = 0.25;
    double e_5_0 = 0.0;
    double e_6_0 = 0.0;
    double e_7_0 = 0.0;
    double e_8_0 = 0.0;
    double e_1_1 = 1.0 / 3.0;
    double e_2_1 = 0.125;
    double e_3_1 = 0.25;
    double e_4_1 = 1.0 / 6.0;
    double e_5_1 = 1.0 / 16.0;
    double e_6_1 = 0.0;
    double e_7_1 = 1.0 / 12.0;
    double e_8_1 = 1.0 / 32.0;
    double e_2_2 = 0.2;
    double e_3_2 = e_0_5;
    double e_4_2 = e_1_5;
    double e_5_2 = 0.1;
    double e_6_2 = 0.0;
    double e_7_2 = 0.0;
    double e_8_2 = 0.0;
    double e_3_3 = 1.0 / 3.0;
    double e_4_3 = 1.0 / 6.0;
    double e_5_3 = 0.0;
    double e_6_3 = 0.125;
    double e_7_3 = 1.0 / 16.0;
    double e_8_3 = e_1_8;
    double e_4_4 = 1.0 / 9.0;
    double e_5_4 = 1.0 / 24.0;
    double e_6_4 = 1.0 / 16.0;
    double e_7_4 = 1.0 / 24.0;
    double e_8_4 = 1.0 / 64.0;
    double e_5_5 = 1.0 / 15.0;
    double e_6_5 = 0.025;
    double e_7_5 = 1.0 / 64.0;
    double e_6_6 = 0.2;
    double e_7_6 = 0.1;
    double e_8_6 = 0.0;
    double e_7_7 = 1.0 / 15.0;
    double e_8_7 = 0.0;
    double e_8_8 = 0.04;



    double x = 0.5;
    double y = 0.5;
    double xy = 0.25;
    double x_2 = 1.0 / 3.0; //x^2
    double y_2 = 1.0 / 3.0; //y^2
    double x_2_y = 1.0 / 6.0; //x^2 * y
    double x_y_2 = 1.0 / 6.0; //x* y^2
    double x_2_y_2 = 1.0 / 9.0; // x^2 * y^2

    Eigen::Matrix4d check_M = Eigen::Matrix4d::Zero();


    auto M_dense = M.makeDense();


    EXPECT_TRUE(M_dense.isApprox(check_M, std::numeric_limits<double>::epsilon()));
    
}

// TEST(dgfe_O1_providers, O1massMatrixAnd01LocalLoadVector){

//     //MESH FROM FILE -----------------------------
//     std::filesystem::path here = __FILE__;
//     auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_100_cells.vtk";
//     lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
//     auto mesh_ptr = reader.mesh();

//     //SIMPLE TEST MESH ----------------------------------
//     //auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

//     //UNIT SQUARE SINGLE POLYGON MESH--------------------
//     // using coord_t = Eigen::Vector2d;
//     // using size_type = lf::mesh::Mesh::size_type;
//     // double scale = 1.0;
//     // std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
//     // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
//     // mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3}}, nullptr);
//     // auto mesh_ptr = mesh_factory_ptr->Build();


//     //setup of dofhandler
//     std::map<lf::base::RefEl, base::size_type> dofmap;
//     dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 9));
//     lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

//     //galerkin matrix initialization
//     lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
//     M.setZero();

//     //initialization of element matrix provider
//     lf::dgfe::DGFEO2MassElementMatrix massMatrixProvider;

//     //assemble mass matrix
//     lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

//     //initialize polynomial for the rhs load vector: 2 + 1.5(x^2 * y)
//     lf::dgfe::Polynomial polynomial;
//     polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
//     polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
//     //initialize load vector provider
//     DGFEO2LocalLoadVector vectorProvider(polynomial);

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
//     auto solution_lambda = [](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
//         return (2.0 + 1.5 * x[0] * x[0] * x[1]);
//     };

//     //setup dgfe space
//     std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
//     //setup mesh function
//     lf::dgfe::MeshFunctionDGFE<double> sol_mesh_function(dgfe_space_ptr, sol_vec);
    
//     //compute l2 error
//     double l2_barycenter_error = lf::dgfe::L2ErrorSubTessellation(sol_mesh_function, solution_lambda, 5);

//     std::cout << "The error of the solution via MeshFunction is: " << l2_barycenter_error << "\n";
// }

// TEST(dgfe_O2_providers, O2massMatrixAnd02LocalLoadVector){

//     //MESH FROM FILE -----------------------------
//     // std::filesystem::path here = __FILE__;
//     // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_100_cells.vtk";
//     // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
//     // auto mesh_ptr = reader.mesh();

//     //SIMPLE TEST MESH ----------------------------------
//     auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

//     //UNIT SQUARE SINGLE POLYGON MESH--------------------
//     // using coord_t = Eigen::Vector2d;
//     // using size_type = lf::mesh::Mesh::size_type;
//     // double scale = 1.0;
//     // std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);
//     // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 0.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 0.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({1.0 * scale, 1.0 * scale}));
//     // mesh_factory_ptr->AddPoint(coord_t({0.0 * scale, 1.0 * scale}));
//     // mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3}}, nullptr);
//     // auto mesh_ptr = mesh_factory_ptr->Build();


//     //setup of dofhandler
//     std::map<lf::base::RefEl, base::size_type> dofmap;
//     dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 9));
//     lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

//     //galerkin matrix initialization
//     lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
//     M.setZero();

//     //initialization of element matrix provider
//     lf::dgfe::DGFEO2MassElementMatrix massMatrixProvider;

//     //assemble mass matrix
//     lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

//     //initialize polynomial for the rhs load vector: 2 + 1.5(x^2 * y)
//     lf::dgfe::Polynomial polynomial;
//     polynomial.push_back(std::make_pair(2.0, std::make_pair(0,0)));
//     polynomial.push_back(std::make_pair(1.5, std::make_pair(2,1)));
//     //initialize load vector provider
//     lf::dgfe::DGFEO2LocalLoadVector vectorProvider(polynomial);

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

//     //setup dgfe space
//     std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
//     //setup mesh function
//     lf::dgfe::MeshFunctionDGFE<double> sol_mesh_function(dgfe_space_ptr, sol_vec);
    
//     //compute l2 error
//     double l2_barycenter_error = lf::dgfe::L2ErrorBarycenter(sol_mesh_function, solution_lambda);

//     std::cout << "The error of the solution via MeshFunction is: " << l2_barycenter_error << "\n";
// }

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