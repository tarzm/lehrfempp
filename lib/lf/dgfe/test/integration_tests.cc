/**
 * @file
 * @brief Tests of the integration algorithm
 * @author Tarzis Maurer
 * @date 2022-04-19
 * @copyright ETH Zurich
 */

#include <cmath>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>

#define NORMALTOLERANCE 0.0000000001

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

    //check outwardNBormal of all edges of the pentagon from the paper
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
    EXPECT_NEAR(euclideanDist(b_point, c_point), 3.2015621187164243, 1e-10);
    Eigen::MatrixXd d_point(2,1);
    Eigen::MatrixXd e_point(2,1);
    d_point << 4.1, -1.7;
    e_point << 0.8, -0.9;
    EXPECT_NEAR(euclideanDist(d_point, e_point), 3.39559, 1e-5);
}

TEST(integration, triangle){

    //test simple triangle from paper
    Eigen::MatrixXd a_polygon(2,3);
    a_polygon <<    -1.0, 1.0, -1.0,
                    -1.0, 0.0, 1.0;
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 10), 0.0111339078, 1e-10);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 5), 0.0, 1e-10);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 20), -0.005890191, 1e-10);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 20, 20), 0.0030396808, 1e-10);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 40, 40), 7.9534562047e-14, 1e-15);
}

TEST(integration, pentagon){

    //test pentagon from paper
    Eigen::MatrixXd b_polygon(2,5);
    b_polygon <<    -0.666666666666667, 0.555555555555556, 1.0, -0.555555555555556, -1.0,
                    -0.789473684210526, -1.0, -0.52631578947368, 1.0, -0.157894736842105;
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 5), -0.0020324991, 1e-8);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 10), 7.4274779926e-5, 1e-8);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 5), -2.0911953867e-4, 1e-8);

    
}

TEST(integration, nonConvexPolygon){
    //test non-convex polygon from paper
    Eigen::MatrixXd c_polygon(2,15);
    c_polygon <<    0.413048522141662, 0.024879797655533, -0.082799691823524, -0.533191422779328, -0.553573605852999, -0.972432940212767, -1.000000000000000, -0.789986179147920, -0.627452906935866, -0.452662174765764, -0.069106265580153, 0.141448047807069, 1.000000000000000, 0.363704451489016, 0.627086024018283,
                    0.781696234443715, 0.415324992429711, 0.688810136531751, 1.000000000000000, 0.580958514816226, 0.734117068746903, 0.238078507228890, 0.012425068086110, -0.636532897516109, -1.000000000000000, -0.289054989277619, -0.464417038155806, -0.245698820584615, -0.134079689960635, -0.110940423607648;
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 5), -0.002589861, 1e-10);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 10, 5), 0.0014996521, 1e-10);


}

} //namespace lf::dgfe::test