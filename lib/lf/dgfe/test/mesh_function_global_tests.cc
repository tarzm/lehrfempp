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


} //namespace lf::dgfe::test