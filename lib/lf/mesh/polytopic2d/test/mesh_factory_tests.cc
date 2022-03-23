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



//checks whether two PolygonPairs have the same entries neglecting the order
bool CheckEqualPolygonPair(PolygonPair polypair1, PolygonPair polypair2){
    return ((polypair1.first == polypair2.first && polypair1.second == polypair2.second) || (polypair1.first == polypair2.second && polypair1.second == polypair2.first));
}

TEST(lf_mesh_factory, EdgePolygonAdjacency){
    //optain test mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //set up Adjacency MeshDataSet
    auto adjacency = lf::mesh::polytopic2d::EdgePolygonAdjacency(mesh_ptr);

    //check some adjacencies by hand
    auto poly0 = mesh_ptr->Entities(0)[0];
    auto poly3 = mesh_ptr->Entities(0)[3];
    auto poly5 = mesh_ptr->Entities(0)[5];

    auto pair_0_3 = std::make_pair(poly0, poly3);
    auto pair_5_3 = std::make_pair(poly5, poly3);
    auto pair_3_5 = std::make_pair(poly3, poly5);
    auto pair_0_null = std::make_pair(poly0, nullptr);
    auto pair_null_0 = std::make_pair(nullptr, poly0);
    auto pair_null_5 = std::make_pair(nullptr, poly5);

    auto edge4 = mesh_ptr->Entities(1)[4];
    auto edge2 = mesh_ptr->Entities(1)[2];
    auto edge12 = mesh_ptr->Entities(1)[12];
    auto edge17 = mesh_ptr->Entities(1)[17];

    //First test CheckEqualPolygonPair function
    EXPECT_TRUE(CheckEqualPolygonPair( pair_5_3 , pair_3_5) );
    EXPECT_FALSE(CheckEqualPolygonPair( pair_0_3 , pair_3_5) );

    //Now test the CodimMeshDataSet EdgePolygonAdjacency
    EXPECT_TRUE(CheckEqualPolygonPair( pair_0_3 , adjacency(*edge2) ));
    EXPECT_TRUE(CheckEqualPolygonPair( pair_0_null , adjacency(*edge4) ));
    EXPECT_TRUE(CheckEqualPolygonPair( pair_5_3 , adjacency(*edge12) ));
    EXPECT_TRUE(CheckEqualPolygonPair( pair_null_5 , adjacency(*edge17) ));
    EXPECT_FALSE(CheckEqualPolygonPair( adjacency(*edge2) , adjacency(*edge12)) );

}

}   //  namespace lf::mesh::polytopic2d