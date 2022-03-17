/**
 * @file
 * @brief Tests for the lf::mesh::polytopic2d::MeshFactory
 * @author Tarzis Maurer
 * @date   2022-02-14
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>

#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

#include <iostream>

#include "mesh_factory_tests.h"

namespace lf::mesh::polytopic2d::test{

TEST(lf_mesh_factory, PolygonConstructor){

    //Coordinates of the Points of the small test mesh
    Eigen::MatrixXd coords(2,6);
    coords <<   0.0, 1.0, 1.0, 0.5, 0.0, 2.0,
                0.0, 0.0, 1.0, 1.5, 1.0, 2.0;

    auto factory = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2,true);
    //Add all Points to the MeshFactory
    for (int i = 0; i < coords.cols(); i++){
        Eigen::Vector2d coords_point = coords.col(i);
        factory->AddPoint(coords_point);
    }

    //Add two Polygons to the MeshFactory, second has the first edge in negative orientation
    factory->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3,4}}, nullptr);
    factory->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,3>{{3,2,5}}, nullptr);

    //Build the mesh
    auto mesh_ptr = factory->Build();

    //Coords(Polygon with index 0) is equal to the coordinates of the Points it has been initialized with
    EXPECT_EQ(lf::mesh::polytopic2d::Corners(mesh_ptr->Entities(0)[0]), coords.block(0, 0, 2, 5));
}

TEST(lf_mesh_factory, IncompleteMesh){
    //Coordinates of the Points of the small test mesh
    Eigen::MatrixXd coords(2,6);
    coords <<   0.0, 1.0, 1.0, 0.5, 0.0, 2.0,
                0.0, 0.0, 1.0, 1.5, 1.0, 2.0;

    auto factory = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2,true);
    //Add all Points to the MeshFactory
    for (int i = 0; i < coords.cols(); i++){
        Eigen::Vector2d coords_point = coords.col(i);
        factory->AddPoint(coords_point);
    }

    //Add only one Polygon to the MeshFactory, node 5 is not used
    factory->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3,4}}, nullptr);

    //Build the mesh
    EXPECT_DEATH(factory->Build(), "Mesh is incomplete");
}

TEST(lf_mesh_factory, Orientations){
    
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //vector which will contain the pairs (polygon_index, local_edge_index) of all edges which have negative orientation within the polygon
    std::vector<std::pair<size_type, size_type>> negative_orientations;

    for(auto &cell : mesh_ptr->Entities(0)){
        auto polygon = dynamic_cast<const lf::mesh::polytopic2d::Polygon*>(cell);
        int counter = 0;
        for (auto &ori : polygon->RelativeOrientations()){
            if (ori == lf::mesh::Orientation::negative){
                negative_orientations.emplace_back(polygon->index(), counter);
            }
            counter ++;
        }
    }

    //Set true solutions by hand
    std::vector<std::pair<size_type, size_type>> negative_orientations_check;
    negative_orientations_check.emplace_back(1,3);
    negative_orientations_check.emplace_back(2,4);
    negative_orientations_check.emplace_back(3,0);
    negative_orientations_check.emplace_back(3,1);
    negative_orientations_check.emplace_back(3,4);
    negative_orientations_check.emplace_back(4,0);
    negative_orientations_check.emplace_back(4,1);
    negative_orientations_check.emplace_back(5,0);
    negative_orientations_check.emplace_back(5,1);
    negative_orientations_check.emplace_back(5,4);

    EXPECT_EQ(negative_orientations, negative_orientations_check);
}

}   //  namespace lf::mesh::polytopic2d