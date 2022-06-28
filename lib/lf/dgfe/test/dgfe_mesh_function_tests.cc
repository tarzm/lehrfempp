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
    auto true_sol_lambda = [](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
        return x[0];
    };

    //setup dgfe space with legendre polynomials of degree 2
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));

    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vector);

    //calculate with mesh function error function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation(dgfe_mesh_function, true_sol_lambda, 2);

    //mesh function should be exact
    EXPECT_NEAR(0.0, mesh_func_l2_error, std::numeric_limits<double>::epsilon());
}

} // namespace lf::dgfe::test