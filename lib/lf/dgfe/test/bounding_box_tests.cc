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




} //namespace lf::dgfe::test