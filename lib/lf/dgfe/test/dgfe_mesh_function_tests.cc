/**
 * @file
 * @brief Tests of the dgfe mesh function functionalities
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
 */

#include <cmath>
#include <filesystem>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include <lf/assemble/assemble.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(meshFunction, O1L2ErrorSubTessellation){
    //retrieve mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(1,1);

    //setup dof vector by hand
    Eigen::VectorXd dof_vector = Eigen::VectorXd::Zero(4*9);
    //cell 0
    dof_vector[0] = 0.25;
    dof_vector[3] = 1;
    //cell 1
    dof_vector[9 + 0] = 0.75;
    dof_vector[9 + 3] = 1;
    //cell 2
    dof_vector[2*9 + 0] = 0.25;
    dof_vector[2*9 + 3] = 1;
    //cell 3
    dof_vector[3*9 + 0] = 0.75;
    dof_vector[3*9 + 3] = 1;

    //setup dgfe space with legendre polynomials of degree 2
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
    //setup mesh function
    lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vector);

    //true solution is 0.25
    auto true_sol_lambda = [](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
        return 0.25;
    };

    //calculate with mesh function
    double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation(dgfe_mesh_function, true_sol_lambda, 9);

    //mesh function should be exact
    EXPECT_EQ(0.0, mesh_func_l2_error);

    for(auto polygon : mesh_ptr->Entities(0)){
        auto box = lf::dgfe::BoundingBox(*polygon);

        Eigen::MatrixXd corner_points(2,4);
        corner_points <<    0,1,1,0,
                            0,0,1,1;


        for (int i = 0; i < 4; i++){
            std::cout << "Polygon " << mesh_ptr->Index(*polygon) << " at point " << i << " has f = " << dgfe_mesh_function(*polygon, corner_points.col(i))[0] << "\n";
        }
    }



}



// TEST(meshFunction, l2BarycenterError){

//     //retrieve mesh
//     auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(1,1);

//     //setup dof vector by hand
//     Eigen::VectorXd dof_vector = Eigen::VectorXd::Zero(4*9);
//     //cell 0
//     dof_vector[4] = 1.5;
//     dof_vector[8] = 0.8;
//     //cell 1
//     dof_vector[9 + 5] = 2.7;
//     //cell 2
//     dof_vector[2*9 + 6] = 1.1;
//     //cell 3
//     dof_vector[3*9 + 7] = -0.5;

//     //setup dgfe space with legendre polynomials of degree 2
//     std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, 2));
//     //setup mesh function
//     lf::dgfe::MeshFunctionDGFE<double> dgfe_mesh_function(dgfe_space_ptr, dof_vector);

//     // dgfe function of cell 0
//     auto dgfe_lambda_0 = [](Eigen::Vector2d x) -> double {
//         std::cout << "Lambda Term 0: " << x[0] * x[1] << " and " << (1.5 * x[0] * x[0] - 0.5) * (1.5 * x[1] * x[1] - 0.5) << "\n";
        
//         return (1.5 * x[0] * x[1] + 0.8 * (1.5 * x[0] * x[0] - 0.5) * (1.5 * x[1] * x[1] - 0.5));
//     };
//     // dgfe function of cell 1
//     auto dgfe_lambda_1 = [](Eigen::Vector2d x) -> double {
//         //std::cout << "Lambda Term 1: " << x[0] * (1.5 * x[1] * x[1] - 0.5) << "\n";
//         std::cout << "Legendre " << multiIndexToDegree(5, 2).first << " is " << legendre_polynomial(multiIndexToDegree(5, 2).first, x[0]);
//         std::cout << ", Legendre " << multiIndexToDegree(5, 2).second << " is " << legendre_polynomial(multiIndexToDegree(5, 2).first, x[1]) << "\n";
//         return (2.7 * x[0] * (1.5 * x[1] * x[1] - 0.5));
//     };
//     // dgfe function of cell 2
//     auto dgfe_lambda_2 = [](Eigen::Vector2d x) -> double {
//         //std::cout << "Lambda Term 2: " << x[0] * (1.5 * x[1] * x[1] - 0.5) << "\n";
//         std::cout << "Legendre " << multiIndexToDegree(6, 2).first << " is " << legendre_polynomial(multiIndexToDegree(6, 2).first, x[0]);
//         std::cout << ", Legendre " << multiIndexToDegree(6, 2).second << " is " << legendre_polynomial(multiIndexToDegree(6, 2).first, x[1]) << "\n";
//         return (1.1 * (1.5 * x[0] * x[0] - 0.5));
//     };
//     // dgfe function of cell 3
//     auto dgfe_lambda_3 = [](Eigen::Vector2d x) -> double {
//         //std::cout << "Lambda Term 3: " << (1.5 * x[0] * x[0] - 0.5) * x[1] << "\n";
//         std::cout << "Legendre " << multiIndexToDegree(7, 2).first << " is " << legendre_polynomial(multiIndexToDegree(7, 2).first, x[0]);
//         std::cout << ", Legendre " << multiIndexToDegree(7, 2).second << " is " << legendre_polynomial(multiIndexToDegree(7, 2).first, x[1]) << "\n";
//         return (-0.5 * (1.5 * x[0] * x[0] - 0.5) * x[1]);
//     };

//     //Coordinates of the mapped barycenters
//     Eigen::MatrixXd mapped_bary_coords(2,1);
//     mapped_bary_coords <<  0.0, 0.0;

//     //true solution is 3 + 1.2 x^2 * y
//     auto true_sol_lambda = [](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
//         return (3.0 + 1.2 * x[0] * x[0] * x[1]);
//     };

//     //compute error by hand
//     double error_check = 0;

//     error_check += std::abs(dgfe_lambda_0(mapped_bary_coords.col(0)) - true_sol_lambda(nullptr, mapped_bary_coords.col(0))) * /* area */ 0.25;
//     // std::cout << "Cell 0 error: " << std::abs(dgfe_lambda_0(mapped_bary_coords.col(0)) - true_sol_lambda(mapped_bary_coords.col(0))) * /* area */ 0.25 << "\n";
//     // std::cout << "Eval " << 0 << ": " << dgfe_lambda_0(mapped_bary_coords.col(0)) << "\n";

//     error_check += std::abs(dgfe_lambda_1(mapped_bary_coords.col(0)) - true_sol_lambda(nullptr, mapped_bary_coords.col(0))) * /* area */ 0.25;
//     // std::cout << "Cell 1 error: " << std::abs(dgfe_lambda_1(mapped_bary_coords.col(1)) - true_sol_lambda(mapped_bary_coords.col(1))) * /* area */ 0.25 << "\n";
//     // std::cout << "Eval " << 1 << ": " << dgfe_lambda_1(mapped_bary_coords.col(1)) << "\n";

//     error_check += std::abs(dgfe_lambda_2(mapped_bary_coords.col(0)) - true_sol_lambda(nullptr, mapped_bary_coords.col(0))) * /* area */ 0.25;
//     // std::cout << "Cell 2 error: " << std::abs(dgfe_lambda_2(mapped_bary_coords.col(2)) - true_sol_lambda(mapped_bary_coords.col(2))) * /* area */ 0.25 << "\n";
//     // std::cout << "Eval " << 2 << ": " << dgfe_lambda_2(mapped_bary_coords.col(2)) << "\n";

//     error_check += std::abs(dgfe_lambda_3(mapped_bary_coords.col(0)) - true_sol_lambda(nullptr, mapped_bary_coords.col(0))) * /* area */ 0.25;
//     // std::cout << "Cell 3 error: " << std::abs(dgfe_lambda_3(mapped_bary_coords.col(3)) - true_sol_lambda(mapped_bary_coords.col(3))) * /* area */ 0.25 << "\n";
//     // std::cout << "Eval " << 3 << ": " << dgfe_lambda_3(mapped_bary_coords.col(3)) << "\n";

//     //calculate with mesh function
//     double mesh_func_l2_error = lf::dgfe::L2ErrorSubTessellation(dgfe_mesh_function, true_sol_lambda, 4);

//     EXPECT_EQ(error_check, mesh_func_l2_error);
// }






} // namespace lf::dgfe::test