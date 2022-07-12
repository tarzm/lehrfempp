/**
 * @file
 * @brief Tests of the discontinuity penalization functionalities
 * @author Tarzis Maurer
 * @date June 22
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

auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(1,1);

for (auto cell : mesh_ptr->Entities(0)){
    for (auto edge : mesh_ptr->)
}



} //namespace lf::dgfe::test