/**
 * @file
 * @brief Tests of the BoundingBox class
 * @author Tarzis Maurer
 * @date 2022-05-03
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
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    for (auto cell : mesh_ptr->Entities(0)){
        BoundingBox box(*cell);
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        auto inverse_mapped = box.inverseMap(corners);
        auto mapped = box.map(inverse_mapped);
        EXPECT_TRUE(corners.isApprox(mapped));
    }
}

TEST(boundingBox, basic2){

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

    auto polygon_ptr = mesh_ptr->Entities(0)[0];

    BoundingBox box(*polygon_ptr);
    auto corners = lf::mesh::polytopic2d::Corners(polygon_ptr);
    auto inverse_mapped = box.inverseMap(corners);

    Eigen::MatrixXd check_map(2,4);
    check_map <<    -1,1,1,-1,
                    -1,-1,1,1;
    
    EXPECT_EQ(inverse_mapped, check_map);
}




} //namespace lf::dgfe::test