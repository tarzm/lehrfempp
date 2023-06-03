/**
 * @file
 * @brief Tests of the dgfe global mesh function functionalities
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

TEST(meshFunctionGlobal, basic){

    //retrieve mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(1,1);

    // mesh function f(x) := e^(x*y^2) * sin(y * PI)    
    auto lambda = [](Eigen::Vector2d x) -> double {
        return exp(x[0] * x[1] * x[1]) * std::sin(x[1] * M_PI);
    };
    lf::dgfe::MeshFunctionGlobalDGFE<decltype(lambda)> mesh_function(lambda);

    Eigen::MatrixXd barycenters(2,4);
    barycenters <<  0.25,   0.75,   0.25,   0.75,
                    0.25,   0.25,   0.75,   0.75;

            
    //lambda evaluated at the barycenter of each cell
    Eigen::Vector4d check_solution(4);
    for (int i = 0; i < 4; i++){
        check_solution[i] = lambda(barycenters.col(i)); 
    }

    Eigen::Vector2d zero_point{0.0, 0.0};

    //check each entry
    for (auto cell : mesh_ptr->Entities(0)){
        EXPECT_NEAR(mesh_function(*cell, zero_point)[0], check_solution[mesh_ptr->Index(*cell)], std::numeric_limits<double>::epsilon());
    }
}

TEST(L2ErrorSubTessellation, twoGlobalMeshFunctions){

    //retrieve mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    auto exp_lambda = [](Eigen::Vector2d x) -> double {
        return exp(x[0]);
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_exp{exp_lambda};

    auto sin_lambda = [](Eigen::Vector2d x) -> double {
        return sin(x[0] * x[1] * M_PI);
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_sin{sin_lambda};

    //calculated via https://www.wolframalpha.com/input?i=integrate%5B+%28+exp%28x%29+-+sin%28x*y*PI%29%29%5E2%2C+%7Bx%2C+0%2C+1%7D%2C+%7By%2C+0%2C+1%7D%5D
    double check = std::sqrt(1.601562432912688047891547320585505904407596611809606829609916392874169123648991989538269064291911348);

    double l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(m_exp), decltype(m_sin)>(m_exp, m_sin, mesh_ptr, 20);

    EXPECT_NEAR( l2_error, check, std::numeric_limits<double>::epsilon());

}

TEST(L2ErrorSubTessellation, DGFEandGlobalMeshFunction){

    //retrieve mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 1);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);
    //dofhandler
    auto dof_handler = dgfe_space_ptr->LocGlobMap();

    //create dof vector by hand => Mesh Funtion is constant 2.5
    auto n_dofs = dof_handler.NumDofs();
    Eigen::VectorXd dof_vec(n_dofs);
    dof_vec.setZero();
    for (int i = 0; i < n_dofs; i += 4){
        dof_vec[i] = 2.5;
    }

    //Mesh_function
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vec);

    auto sin_lambda = [](Eigen::Vector2d x) -> double {
        return sin(x[0] * x[1] * M_PI);
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_sin{sin_lambda};

    //calculated via https://www.wolframalpha.com/input?i=integrate%5B+%28+2.5+-+sin%28x*y*PI%29%29%5E2%2C+%7Bx%2C+0%2C+1%7D%2C+%7By%2C+0%2C+1%7D%5D
    double check = 4.013831745425869723500183307733026855967725529923528045789347059812650762131396931799738047219165127;

    double l2_error = lf::dgfe::L2ErrorSubTessellation<double, decltype(dgfe_mesh_function), decltype(m_sin)>(dgfe_mesh_function, m_sin, mesh_ptr, 20);

    EXPECT_NEAR( l2_error, check, std::numeric_limits<double>::epsilon());

}

TEST(L2ErrorSubTessellation, twoGlobalMeshFunctionsGrad){

    //retrieve mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    auto u_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d() << std::sin(x[0]), std::sin(x[1])).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_u{u_lambda};

    auto g_lambda = [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d() << exp(x[1]), exp(x[0])).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_g{g_lambda};

    //calculated via https://www.wolframalpha.com/input?i=integrate%5B+%28+sin%28x%29+-+exp%28y%29%29%5E2+%2B+%28sin%28y%29+-+exp%28x%29%29%5E2%2C+%7Bx%2C+0%2C+1%7D%2C+%7By%2C+0%2C+1%7D%5D
    double check = std::sqrt(3.774846607872610197043391601742996031178784974072409870540084489580129908234637109607810992350851787);

    EXPECT_NEAR( lf::dgfe::L2ErrorGradSubTessellation<double>(m_u, m_g, mesh_ptr, 20), check, 10 * std::numeric_limits<double>::epsilon());

}

} //namespace lf::dgfe::test