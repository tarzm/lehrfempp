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

TEST(discontinuity_penalization, simplexArea){
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(1,1);

    for (auto cell : mesh_ptr->Entities(0)){
        unsigned edge_loc_idx = 0;

        for (auto edge : cell->SubEntities(1)){
            auto simplex_areas = lf::dgfe::simplexAreas(*cell, edge_loc_idx);

            for (int i = 0; i < simplex_areas.size(); i++){
                EXPECT_EQ(simplex_areas[i], 0.125);
            }
            edge_loc_idx++;
        }
    }

}

TEST(discontinuity_penalization, basicEvaluation){
    //mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    lf::dgfe::DiscontinuityPenalization disc_pen(dgfe_space_ptr, 1.0, 1.0);
    double A_f = 1.0;

    for (auto edge : mesh_ptr->Entities(1)){
        double alpha = disc_pen(*edge, A_f);
    }
}


} //namespace lf::dgfe::test