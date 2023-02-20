/**
 * @file
 * @brief Tests of the sub tessellation functionalities
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include <lf/assemble/assemble.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {


TEST(sub_tessellation, basic){

    //UNIT SQUARE SINGLE POLYGON MESH, only one polygon which is the square itself
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

    auto polygon = mesh_ptr->Entities(0)[0];

    auto sub_tessellation = subTessellation(polygon);
    auto corners = lf::mesh::polytopic2d::Corners(polygon);
    Eigen::Vector2d barycenter;
    barycenter << 0.5, 0.5;

    int node_counter = 0;
    //loop over created triangle geometries
    auto tessellation_size = sub_tessellation.size();
    for (int i = 0; i < tessellation_size; i++){
        Eigen::MatrixXd check_tria_corners(2,3);
        check_tria_corners.col(0) = corners.col(node_counter);
        check_tria_corners.col(1) = corners.col((node_counter + 1) % 4);
        check_tria_corners.col(2) = barycenter;
        auto tria_corners = lf::geometry::Corners(*sub_tessellation[i]);
        EXPECT_EQ(tria_corners, check_tria_corners);
        node_counter++;
    }
}






} //namespace lf::dgfe::test