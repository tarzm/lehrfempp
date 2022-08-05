/**
 * @file
 * @brief Tests of the auxiliary operators functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
 */


#include <gtest/gtest.h>

#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(boundingBox, basic){

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

    auto cell = *(mesh_ptr->Entities(0)[0]);

    // 2x2 diffusion tensor A(x)
    auto a_simple_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
        return (Eigen::Matrix<double, 2, 2>() << 1.0, 0.0, 0.0, 1.0).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_simple{a_simple_lambda};

    //sub-tessellation
    auto subtessellation = lf::dgfe::subTessellation(cell);
    //box
    lf::dgfe::BoundingBox box(cell);

    //quadrule setup
    auto qr_t = make_quadrule(lf::base::RefEl::kTria(), 10);
    // qr points
    const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
    //weights
    Eigen::VectorXd w_ref_t{qr_t.Weights()};

    




}

}