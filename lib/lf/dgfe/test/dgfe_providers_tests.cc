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
#include <typeinfo>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include <lf/assemble/assemble.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(dgfe_subTessellation_providers, singleSquaremassMatrixProviderO1){
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
    lf::dgfe::DGFEMassElementMatrixST<double> massMatrixProvider(7, 1); //quadrule of degree 4 would be ok too but then coefficients inside the matrix would be bigger than epsilon
    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //entries of the check matrix => functions integrated over (-1,1)^2 and then divided by 4 (reference box is just 4 times as big as single polygon)
    double x = 0.0;
    double y = 0.0;
    double xy = 0.0;
    double x_2 = 1.0 / 3.0; //x^2
    double y_2 = 1.0 / 3.0; //y^2
    double x_2_y = 0.0; //x^2 * y
    double x_y_2 = 0.0; //x* y^2
    double x_2_y_2 = 1.0 / 9.0; // x^2 * y^2

    Eigen::Matrix4d check_M;
    check_M <<  1.0,      y,      x,      xy,
                y,      y_2,    xy,     x_y_2,
                x,      xy,     x_2,    x_2_y,
                xy,     x_y_2,  x_2_y,  x_2_y_2;

    auto M_dense = M.makeDense();

    //compare each entry
    for (int i = 0; i < M_dense.cols(); i++){
        for (int j = 0; j < M_dense.cols(); j++){
            EXPECT_NEAR(M_dense(i,j), check_M(i,j), std::numeric_limits<double>::epsilon()) << "At (" << i << ", " << j << ")\n";
        }
    }
}

TEST(dgfe_subTessellation_providers, singleSquaremassMatrixProviderO2){

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
    lf::dgfe::DGFEMassElementMatrixST<double> massMatrixProvider(9, 2); 
    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    /**entries of the check matrix => functions integrated over (0,1)^2
     *BASIS FUNCTION IDX  | FUNCTION EXPRESSION
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
     * Now the entry (i, j) of the elements matrix is basis_i(x,y) * basis_j(x,y) integrated over (-1,1)^2 and divided by 4
    */
   
    double e_0_0 = 1.0;
    double e_1_0 = 0.0;
    double e_2_0 = 0.0;
    double e_3_0 = 0.0;
    double e_4_0 = 0.0;
    double e_5_0 = 0.0;
    double e_6_0 = 0.0;
    double e_7_0 = 0.0;
    double e_8_0 = 0.0;

    double e_1_1 = 1.0 / 3.0;
    double e_2_1 = 0.0;
    double e_3_1 = 0.0;
    double e_4_1 = 0.0;
    double e_5_1 = 0.0;
    double e_6_1 = 0.0;
    double e_7_1 = 0.0;
    double e_8_1 = 0.0;

    double e_2_2 = 0.2;
    double e_3_2 = 0.0;
    double e_4_2 = 0.0;
    double e_5_2 = 0.0;
    double e_6_2 = 0.0;
    double e_7_2 = 0.0;
    double e_8_2 = 0.0;

    double e_3_3 = 1.0 / 3.0;
    double e_4_3 = 0.0;
    double e_5_3 = 0.0;
    double e_6_3 = 0.0;
    double e_7_3 = 0.0;
    double e_8_3 = 0.0;

    double e_4_4 = 1.0 / 9.0;
    double e_5_4 = 0.0;
    double e_6_4 = 0.0;
    double e_7_4 = 0.0;
    double e_8_4 = 0.0;

    double e_5_5 = (2.0 / 3.0) / 10.0;
    double e_6_5 = 0.0;
    double e_7_5 = 0.0;
    double e_8_5 = 0.0;

    double e_6_6 = 0.2;
    double e_7_6 = 0.0;
    double e_8_6 = 0.0;

    double e_7_7 = (2.0 / 3.0) / 10.0;
    double e_8_7 = 0.0;

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
    auto load_lambda = [](Eigen::Vector2d x) -> double {
        return exp(x[0] * x[1]);
    };
    lf::dgfe::MeshFunctionGlobalDGFE<decltype(load_lambda)> lambda_msh_funct(load_lambda);

    //rhs vector initialization
    Eigen::VectorXd rhs(dofhandler.NumDofs());
    rhs.setZero();
    //initialization of vector provider
    lf::dgfe::DGFELoadElementVectorProvider<double, decltype(lambda_msh_funct)> vec_provider(dgfe_space_ptr, lambda_msh_funct);
    //assemble load vector
    lf::assemble::AssembleVectorLocally(0, dofhandler, vec_provider, rhs);

    Eigen::VectorXd rhs_check(dofhandler.NumDofs());
    
    //all calculated via wolframalpha.com
    // => integrate {legendreP[n, (x-0.5)/0.5)] * legendreP[n, (y-0.5)/0.5)] * e^(x*y) dx dy} over (0,1) ^ 2
    rhs_check[0] = 1.3179021514544038948600088442492318379749012457928;
    rhs_check[1] = 0.11866150546368657586056609845609315753959294160714;
    rhs_check[2] = 0.0082111807001324826982840161332568514314186835930265;
    rhs_check[3] = 0.11866150546368657586056609845609315753959294160714;
    rhs_check[4] = 0.046293545636788315419973467252304486075296262621648;
    rhs_check[5] = 0.0044144106537190776954977450989502536079602045212707;
    rhs_check[6] = 0.0082111807001324826982840161332568514314186835930265;
    rhs_check[7] = 0.0044144106537190776954977450989502536079602045212707;
    rhs_check[8] = 0.0011434868300297584672265496854779512452763431320375;

    //each coefficient of the vector should be equal up to machine precision
    for (int i = 0; i < 9; i++){
        EXPECT_NEAR(rhs[i], rhs_check(i), std::numeric_limits<double>::epsilon() * 10) << "At " << i << "\n";
    }
}

TEST(dgfe_subTessellation_providers, massMatrixmeshFunction){
    //SIMPLE TEST MESH ----------------------------------
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    // //get mesh
    // std::filesystem::path here = __FILE__;
    // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_400_cells.vtk";
    // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    // auto mesh_ptr = reader.mesh();


    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), 9));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_ptr, dofmap);

    //setup dgfe space
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(dofhandler.NumDofs(), dofhandler.NumDofs());
    M.setZero();

    //initialization of element matrix provider
    lf::dgfe::DGFEMassElementMatrixST<double> massMatrixProvider(10, 2);

    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dofhandler, dofhandler, massMatrixProvider, M);

    //assemble rhs vector
    //load function is e^(x*y)
    auto load_lambda = [](Eigen::Vector2d x) -> double {
        return exp(x[0] * x[1]);
    };
    lf::dgfe::MeshFunctionGlobalDGFE<decltype(load_lambda)> lambda_msh_funct(load_lambda);

    //rhs vector initialization
    Eigen::VectorXd rhs(dofhandler.NumDofs());
    rhs.setZero();
    //initialization of vector provider
    lf::dgfe::DGFELoadElementVectorProvider<double, decltype(lambda_msh_funct)> vec_provider(dgfe_space_ptr, lambda_msh_funct);
    //assemble load vector
    lf::assemble::AssembleVectorLocally(0, dofhandler, vec_provider, rhs);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(rhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    //setup mesh function with solution vector
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, sol_vec);

    //calculate with mesh function error function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(dgfe_mesh_function), decltype(lambda_msh_funct)>(dgfe_mesh_function, lambda_msh_funct, mesh_ptr, 30);
    
    std::cout << "The error of the solution via MeshFunction is: " << mesh_func_l2_error << "\n";
}

TEST(L2ProjectionSqrtANablaBasisLoadVector, simpleSquare){
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

    //setup dgfe space
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
    unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();

    // 2x2 diffusion tensor A(x)
    auto a_coeff_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
        return (Eigen::Matrix<double, 2, 2>() << 0.0, 0.0, 0.0, 0.0).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{a_coeff_lambda};

    //l2 projection load vector provider setup
    lf::dgfe::L2ProjectionSqrtANablaBasisLoadVector<double, decltype(m_a_coeff)> l2_projection_provider_0(dgfe_space_ptr, m_a_coeff, 0);
    lf::dgfe::L2ProjectionSqrtANablaBasisLoadVector<double, decltype(m_a_coeff)> l2_projection_provider_1(dgfe_space_ptr, m_a_coeff, 1);

    //initialization of element matrix provider
    lf::dgfe::DGFEMassElementMatrixST<double> massMatrixProvider(10, 2); 

    //galerkin matrix initialization
    lf::assemble::COOMatrix<double> M(n_dofs, n_dofs);
    M.setZero();

    //assemble mass matrix
    lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), massMatrixProvider, M);

    //rhs setup
    Eigen::VectorXd rhs_0(n_dofs);
    Eigen::VectorXd rhs_1(n_dofs);
    rhs_0.setZero();
    rhs_1.setZero(); 

    //assemble both rhs vectors
    lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), l2_projection_provider_0, rhs_0);
    lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), l2_projection_provider_1, rhs_1);

    //solve LSE
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec_0 = solver.solve(rhs_0);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
    Eigen::VectorXd sol_vec_1 = solver.solve(rhs_1);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");





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



}