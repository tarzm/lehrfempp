/**
 * @file
 * @brief Tests for the lf::mesh::polytopic2d::Polygon
 * @author Tarzis Maurer
 * @date   2022-02-14
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include "polygon_test.h"

namespace lf::mesh::polytopic2d::test{

TEST(lf_polygon, constructor){
    using vec = Eigen::Vector2d;

    //Filling in some points of the MeshFactory
    Eigen::MatrixXd coords(2,5);
    coords <<   0.0, 1.0, 1.0, 0.5, 0.0,
                0.0, 0.0, 1.0, 1.5, 1.0;

    auto factory = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2,true);
    for (int i = 0; i < 5; i++){
        Eigen::Vector2d coords_point = coords.col(i);
        factory->AddPoint(coords_point);
    }
    
    








    //Constructing a pentagon by hand

    //construct Points
    // std::vector<std::unique_ptr<lf::geometry::Point>> geometries;
    // std::vector<lf::mesh::hybrid2d::Point*> nodes;

    // Eigen::MatrixXd coords(2,5);
    // coords <<   0.0, 1.0, 1.0, 0.5, 0.0,
    //             0.0, 0.0, 1.0, 1.5, 1.0;

    // //fill nodes vector
    // for (int column = 0; column < 5; column++){
    //     auto point_geo = std::make_unique<geometry::Point>(coords.col(column));
    //     lf::mesh::hybrid2d::Point point(column, std::move(point_geo));
    //     nodes.emplace_back(&point);
    // }

    // for (int column = 0; column < 5; column++){
    //     std::cout << "Hello 0" << "\n";
    //     lf::geometry::Point geo_point= *(nodes.at(column))->Geometry();
    //     std::cout << "Hello 1" << "\n";
    //     // auto node_coord = lf::geometry::Corners(*(node->Geometry()));
    //     // EXPECT_EQ(coords.col(column), node_coord);
    // }
}








}   //  namespace lf::mesh::polytopic2d