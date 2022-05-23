/**
 * @file
 * @brief Test the NumpyPolytopicWriter class
 * @author Tarzis Maurer
 * @date   2022-05-13
 * @copyright ETH Zurich
 */

#include <gtest/gtest.h>
#include <filesystem>

#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::io::test {

TEST(lf_io, numpyPolytopicWriter){

    //retrieve test mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //get mesh
    // std::filesystem::path here = __FILE__;
    // auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_1000_cells.vtk";

    // lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    
    // auto mesh_ptr = reader.mesh();

    // define mesh function of the form sin(\pi*x)*cos(\pi*y) (in global coordinates)
    auto mfTrig = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return std::sin(base::kPi * x[0]) * std::cos(base::kPi * x[1]);
    });

    auto mfSimple = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return x;
    });

    lf::io::NumpyPolytopicWriter writer(mesh_ptr);

    // Eigen::Matrix<double, 0, 1> dummy;

    // auto node_0 = mesh_ptr->Entities(2)[4];

    // std::cout << "Value of meshFunction node 4 is " << mfTrig(*node_0, dummy)[0] << "\n";

    std::cout << "OK 0 OUT\n";

    writer.writeSimple(mfTrig, "testian.txt");



}




} // namespace lf::io::test