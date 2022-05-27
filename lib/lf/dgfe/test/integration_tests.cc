/**
 * @file
 * @brief Tests of the integration algorithm
 * @author Tarzis Maurer
 * @date 2022-04-19
 * @copyright ETH Zurich
 */

#include <cmath>
#include <filesystem>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include "lf/mesh/test_utils/test_meshes.h"

#define NORMALTOLERANCE 1e-12
#define TOLERANCE 1e-10

namespace lf::dgfe::test {

//edge1
bool isNormalOf(const Eigen::MatrixXd normal, const Eigen::MatrixXd edge){
    LF_VERIFY_MSG(std::abs(normal.col(0).norm() - 1.0) < NORMALTOLERANCE, "Normal does not have length 1");
    Eigen::MatrixXd edge_vector(2,1);
    edge_vector(0,0) = edge(0,1) - edge(0,0);
    edge_vector(1,0) = edge(1,1) - edge(1,0);
    return (std::abs(edge_vector.col(0).dot(normal.col(0))) < NORMALTOLERANCE);
}

TEST(integration, helperFunctions){

    //test outwardNormal
    Eigen::MatrixXd a_edge(2,2);
    a_edge <<       0, 1,
                    1, 0;
    Eigen::MatrixXd normal_check(2,1);
    normal_check << - std::sqrt(2.0) / 2.0,  - std::sqrt(2.0) / 2.0;
    EXPECT_TRUE(lf::dgfe::outwardNormal(a_edge).isApprox(normal_check));
    EXPECT_TRUE(isNormalOf(lf::dgfe::outwardNormal(a_edge), a_edge));

    //check outwardNormal of all edges of the pentagon from the paper
    Eigen::MatrixXd b_polygon(2,5);
    b_polygon <<    -0.666666666666667, 0.555555555555556, 1.0, -0.555555555555556, -1.0,
                    -0.789473684210526, -1.0, -0.52631578947368, 1.0, -0.157894736842105;

    for (int i = 0; i < b_polygon.cols(); i++){
        Eigen::MatrixXd edge_i(2,2);
        edge_i.col(0) = b_polygon.col(i);
        edge_i.col(1) = b_polygon.col((i + 1) % b_polygon.cols());
        EXPECT_TRUE(isNormalOf(lf::dgfe::outwardNormal(edge_i), edge_i));
    }

    //test euclideanDist
    Eigen::MatrixXd b_point(2,1);
    Eigen::MatrixXd c_point(2,1);
    b_point << 0.0, 0.0;
    c_point << 2.0, 2.5;
    EXPECT_NEAR(euclideanDist(b_point, c_point), 3.2015621187164243, 1e-12);
    Eigen::MatrixXd f_point(2,1);
    Eigen::MatrixXd g_point(2,1);
    f_point << -0.8, 3.0;
    g_point << 2.34, 1.01;
    EXPECT_NEAR(euclideanDist(f_point, g_point), std::sqrt(std::pow( 2.34-(-0.8) , 2) + std::pow( 1.01-3.0 , 2)), 1e-12);

}

TEST(integration, lineIntegrals){

    Eigen::MatrixXd a_polygon(2,2);
    a_polygon <<    2.0, 5.0,
                    1.0, 3.0;
    EXPECT_DOUBLE_EQ(lf::dgfe::integrate(a_polygon, 4, 3), std::sqrt(13) * (648.0/8.0 + 2700.0/7.0 + 4806.0/6.0 + 4737.0/5.0 + 2792.0/4.0 + 984.0/3.0 + 192.0/2.0 + 16.0));
    
}

TEST(integration, triangle){

    //test simple triangle from paper
    Eigen::MatrixXd a_polygon(2,3);
    a_polygon <<    -1.0, 1.0, -1.0,
                    -1.0, 0.0, 1.0;
    
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 5), 0.0, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 10), 0.0111339078, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 20, 20), 0.0030396808, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 40, 40), 7.9534562047e-14, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 5), 0.0, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 20, 40), 0.0, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 40, 5), 0.0, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 20), -0.005890191, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 40), -0.001868889, TOLERANCE);
}
// 
TEST(integration, pentagon){

    //test pentagon from paper
    Eigen::MatrixXd b_polygon(2,5);
    b_polygon <<    -0.666666666666667, 0.555555555555556, 1.000000000000000, -0.555555555555556, -1.000000000000000,
                    -0.789473684210526, -1.000000000000000, -0.052631578947368, 1.000000000000000, -0.157894736842105;

    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 5), -0.0020324991, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 10), 7.4274779926e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 20, 20), 6.0738145408e-8, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 40, 40), 2.2238524572e-12, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 5), -2.0911953867e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 20, 5), -1.3797380205e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 40, 5), -7.9203571311e-7, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 20), 8.08469022058e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 40), 4.37593748009e-5, TOLERANCE);
}

TEST(integration, nonConvexPolygon){
    //test non-convex polygon from paper
    Eigen::MatrixXd c_polygon(2,15);
    c_polygon <<    0.413048522141662, 0.024879797655533, -0.082799691823524, -0.533191422779328, -0.553573605852999, -0.972432940212767, -1.000000000000000, -0.789986179147920, -0.627452906935866, -0.452662174765764, -0.069106265580153, 0.141448047807069, 1.000000000000000, 0.363704451489016, 0.627086024018283,
                    0.781696234443715, 0.415324992429711, 0.688810136531751, 1.000000000000000, 0.580958514816226, 0.734117068746903, 0.238078507228890, 0.012425068086110, -0.636532897516109, -1.000000000000000, -0.289054989277619, -0.464417038155806, -0.245698820584615, -0.134079689960635, -0.110940423607648;
    //EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 5), -0.002589861, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 10, 10), 1.5738050178e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 20, 20), 1.3793481020e-6, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 40, 40), 4.2588831784e-10, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 10, 5), 0.0014996521, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 20, 5), 7.0356275077e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 40, 5), 2.5065856538e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 20), -1.330384913e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 40), -3.963064075e-5, TOLERANCE);
}

TEST(integration, polytopicTestMesh){
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    lf::dgfe::scalar_t sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        sum += integrate(corners, 3, 4);
    }
    EXPECT_NEAR(sum, 0.05, 1e-14);

    sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        sum += integrate(corners, 5, 7);
    }
    EXPECT_NEAR(sum, (1.0/48.0), 1e-14);

}

// TEST(integration, bigMesh){
//     //get mesh
//     std::filesystem::path here = __FILE__;
//     auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_polytopic_1000_cells.vtk";

//     lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    
//     auto mesh_ptr = reader.mesh();

//     lf::dgfe::scalar_t sum = 0.0;
//     for (auto cell : mesh_ptr->Entities(0)){
//         auto corners = lf::mesh::polytopic2d::Corners(cell);
//         //std::cout << "Cell " << mesh_ptr->Index(*cell) << " contributes " << integrate(corners, 3, 4) << "\n";
//         sum += integrate(corners, 3, 4);
//     }
//     //this mesh has some very small edges => I believe that is where the error is coming from
//     EXPECT_NEAR(sum, 0.05, TOLERANCE);

// }

} //namespace lf::dgfe::test