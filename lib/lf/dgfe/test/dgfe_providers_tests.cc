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

TEST(dgfe_SubTessellation_providers, singleSquaremassMatrixProviderO1){
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

TEST(dgfe_SubTessellation_providers, singleSquaremassMatrixProviderO2){

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
    lf::dgfe::DGFEMassElementMatrixST massMatrixProvider(9, 2); 
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
    double e_7_1 = 0.0;
    double e_8_1 = 0.0;

    double e_2_2 = 0.2;
    double e_3_2 = 0.0;
    double e_4_2 = 1.0 / 16.0;
    double e_5_2 = 0.1;
    double e_6_2 = 0.0;
    double e_7_2 = 0.0;
    double e_8_2 = 0.0;

    double e_3_3 = 1.0 / 3.0;
    double e_4_3 = 1.0 / 6.0;
    double e_5_3 = 0.0;
    double e_6_3 = 0.125;
    double e_7_3 = 1.0 / 16.0;
    double e_8_3 = 0.0;

    double e_4_4 = 1.0 / 9.0;
    double e_5_4 = 1.0 / 24.0;
    double e_6_4 = 1.0 / 16.0;
    double e_7_4 = 1.0 / 24.0;
    double e_8_4 = 1.0 / 64.0;

    double e_5_5 = 1.0 / 15.0;
    double e_6_5 = 0.0;
    double e_7_5 = 1.0 / 64.0;
    double e_8_5 = 0.025;

    double e_6_6 = 0.2;
    double e_7_6 = 0.1;
    double e_8_6 = 0.0;

    double e_7_7 = 1.0 / 15.0;
    double e_8_7 = 0.025;

    double e_8_8 = 0.04;

    Eigen::MatrixXd check_M = Eigen::MatrixXd::Zero(9,9);
    check_M(0,0) = e_0_0;
    check_M(1,0) = e_1_0;
    check_M(2,0) = e_2_0;
    check_M(3,0) = e_3_0;
    check_M(4,0) = e_4_0;
    check_M(5,0) = e_5_0;
    check_M(6,0) = e_6_0;
    check_M(7,0) = e_7_0;
    check_M(8,0) = e_8_0;

    check_M(1,1) = e_1_1;
    check_M(2,1) = e_2_1;
    check_M(3,1) = e_3_1;
    check_M(4,1) = e_4_1;
    check_M(5,1) = e_5_1;
    check_M(6,1) = e_6_1;
    check_M(7,1) = e_7_1;
    check_M(8,1) = e_8_1;

    check_M(2,2) = e_2_2;
    check_M(3,2) = e_3_2;
    check_M(4,2) = e_4_2;
    check_M(5,2) = e_5_2;
    check_M(6,2) = e_6_2;
    check_M(7,2) = e_7_2;
    check_M(8,2) = e_8_2;

    check_M(3,3) = e_3_3;
    check_M(4,3) = e_4_3;
    check_M(5,3) = e_5_3;
    check_M(6,3) = e_6_3;
    check_M(7,3) = e_7_3;
    check_M(8,3) = e_8_3;

    check_M(4,4) = e_4_4;
    check_M(5,4) = e_5_4;
    check_M(6,4) = e_6_4;
    check_M(7,4) = e_7_4;
    check_M(8,4) = e_8_4;

    check_M(5,5) = e_5_5;
    check_M(6,5) = e_6_5;
    check_M(7,5) = e_7_5;
    check_M(8,5) = e_8_5;

    check_M(6,6) = e_6_6;
    check_M(7,6) = e_7_6;
    check_M(8,6) = e_8_6;

    check_M(7,7) = e_7_7;
    check_M(8,7) = e_8_7;

    check_M(8,8) = e_8_8;

    //M is symmetric
    for (int i = 0; i < 8; i++){ //last row is already full
        for (int j = i + 1; j < 9; j++){ //diagonal is full
            check_M(i, j) = check_M(j, i);
        }
    }

    auto M_dense = M.makeDense();

    //compare each entry
    for (int i = 0; i < M_dense.cols(); i++){
        for (int j = 0; j < M_dense.cols(); j++){
            EXPECT_NEAR(M_dense(i,j), check_M(i,j), std::numeric_limits<double>::epsilon()) << "At (" << i << ", " << j << ")\n";
        }
    }
}

TEST(dgfe_subTessellation_providers, singleSquareO2LoadElementVectorST){

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

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    //load function is e^(x*y)
    auto load_lambda = [](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
        return exp(x[0] * x[1]);
    };

    //rhs vector initialization
    Eigen::VectorXd rhs(dofhandler.NumDofs());
    rhs.setZero();
    //initialization of vector provider
    lf::dgfe::DGFELoadElementVectorProvider<double, decltype(load_lambda)> vec_provider(dgfe_space_ptr, load_lambda);
    //assemble load vector
    lf::assemble::AssembleVectorLocally(0, dofhandler, vec_provider, rhs);

    Eigen::VectorXd rhs_check(dofhandler.NumDofs());
    
    //all calculated via wolframalpha.com
    rhs_check[0] = 1.317902151454403894860008844249231837974901245792783992840461196997646107756139482611953646834392207;
    rhs_check[1] = 0.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274;
    rhs_check[2] = 0.09104892427279805256999557787538408101254937710360800357976940150117694612193025869402317658280389627;
    rhs_check[3] = 0.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274;
    rhs_check[4] = 0.4003796770046413405002786271034306597823458479071755821265064307264305225974081119594285316907742200;
    rhs_check[5] = 0.06343634308190952927942505729467500448550581260008085006606474455184673929290481085723564294966714515;
    rhs_check[6] = 0.09104892427279805256999557787538408101254937710360800357976940150117694612193025869402317658280389627;
    rhs_check[7] = 0.06343634308190952927942505729467500448550581260008085006606474455184673929290481085723564294966714515;
    rhs_check[8] = 0.02776699134271494146374838909686999047316899586590587864083636098047417525069836683556001910011406200;

    //each coefficient of the vector should be equal up to machine precision
    for (int i = 0; i < 9; i++){
        EXPECT_NEAR(rhs[i], rhs_check(i), std::numeric_limits<double>::epsilon()) << "\n";
    }
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