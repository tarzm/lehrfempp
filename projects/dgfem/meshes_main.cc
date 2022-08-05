#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <lf/mesh/mesh.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/base/base.h>
#include <lf/dgfe/dgfe.h>


int main(int /*argc*/, char ** /*argv*/) {

    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/meshes/unit_square_polytopic_100_cells.vtk";

    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    
    auto mesh = reader.mesh();

    auto mfTrig = lf::mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return std::sin(lf::base::kPi * x[0]) * std::cos(lf::base::kPi * x[1]);
    });
      
    // output the mesh + mesh function using 1st order cells:
    lf::io::VtkWriter vtk1(mesh, "1storder.vtk");
    vtk1.WritePointData("trig", mfTrig);
    
    // // output the mesh + mesh function using 2nd order cells:
    // VtkWriter vtk2(mesh, "2ndorder.vtk", 0, 2);
    // vtk2.WritePointData("trig", mfTrig);
















































    // std::cout << "Number of codim 0 entities: " << mesh->Entities(0).size() << "\n";

    // //write mesh for matlplotlib
    // lf::io::writeMatplotlib(*mesh, "100.csv");

    // for (auto polygon : mesh->Entities(0)){
    //     auto index = mesh->Index(*polygon);
    //     if (index == 0 || index == 11 || index == 79 || index == 91 || index == 31){
    //         std::cout << "Polygon " << index << " has points ";
    //         for (auto point : polygon->SubEntities(2)){
    //             std::cout << mesh->Index(*point) << " ";
    //         }
    //         std::cout << "\n";
    //     }
    // }

    // for (auto edge : mesh->Entities(1)){
    //     auto index = mesh->Index(*edge);
    //     if (index == 168 || index == 166){
    //         std::cout << "Edge " << index << " has points ";
    //         for (auto point : edge->SubEntities(1)){
    //             std::cout << mesh->Index(*point) << " ";
    //         }
    //         std::cout << "\n";
    //     }
    // }

    // for (auto point : mesh->Entities(2)){
    //     auto index = mesh->Index(*point);
    //     if (index == 27 || index == 69){
    //         std::cout << "Point " << index << " has coordinates" << lf::geometry::Corners(*(point->Geometry())) << "\n";
    //     }
    // }

    return 0;
}