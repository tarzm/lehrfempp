/**
 * @file
 * @brief Test the VtkPolytopicReader
 * @author Tarzis Maurer
 * @date   2022-04-14
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::io::test {

using size_type = mesh::Mesh::size_type;

TEST(lf_io, vtkPolytopicReader100){

    VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2),
                        test_utils::getMeshPath("unit_square_polytopic_100_cells.vtk"));
    
    auto mesh = reader.mesh();

    //test number of cells
    int n_cells = mesh->Entities(0).size();
    EXPECT_EQ(n_cells, 100);

    //test some cells' node indices at random
    std::vector<size_type> cell_57_indices;
    std::vector<size_type> cell_68_indices;

    for (auto node : mesh->Entities(0)[57]->SubEntities(2)){
        cell_57_indices.push_back(mesh->Index(*node));
    }
    for (auto node : mesh->Entities(0)[68]->SubEntities(2)){
        cell_68_indices.push_back(mesh->Index(*node));
    }

    //extracted from python jupyter notebook script by hand
    std::vector<size_type> cell_57_indices_verify{108, 113, 112, 111, 109, 110};
    std::vector<size_type> cell_68_indices_verify{187, 188, 179, 142, 141, 140, 144};

    EXPECT_EQ(cell_57_indices, cell_57_indices_verify);
    EXPECT_EQ(cell_68_indices, cell_68_indices_verify);
}

TEST(lf_io, vtkPolytopicReader1000){

    VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2),
                        test_utils::getMeshPath("unit_square_polytopic_1000_cells.vtk"));
    
    auto mesh = reader.mesh();

    //test number of cells
    int n_cells = mesh->Entities(0).size();
    EXPECT_EQ(n_cells, 1000);

    //test some cells' node indices at random
    std::vector<size_type> cell_57_indices;
    std::vector<size_type> cell_68_indices;
    std::vector<size_type> cell_438_indices;

    for (auto node : mesh->Entities(0)[57]->SubEntities(2)){
        cell_57_indices.push_back(mesh->Index(*node));
    }
    for (auto node : mesh->Entities(0)[68]->SubEntities(2)){
        cell_68_indices.push_back(mesh->Index(*node));
    }
    for (auto node : mesh->Entities(0)[438]->SubEntities(2)){
        cell_438_indices.push_back(mesh->Index(*node));
    }

    //extracted from python jupyter notebook script by hand
    std::vector<size_type> cell_57_indices_verify{635, 77, 1186, 1189, 83, 634};
    std::vector<size_type> cell_68_indices_verify{1500, 1499, 1816, 1813, 1503, 1504};
    std::vector<size_type> cell_438_indices_verify{1896, 1895, 696, 267, 1203, 1202};

    EXPECT_EQ(cell_57_indices, cell_57_indices_verify);
    EXPECT_EQ(cell_68_indices, cell_68_indices_verify);
    EXPECT_EQ(cell_438_indices, cell_438_indices_verify);
}

} //namespace lf::io::test