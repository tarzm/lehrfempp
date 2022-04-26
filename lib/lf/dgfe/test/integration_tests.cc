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

namespace lf::dgfe::test {


TEST(integration, helperFunctions){

    //test outwardNormal
    Eigen::MatrixXd a_edge(2,2);
    a_edge <<       0, 1,
                    1, 0;
    Eigen::MatrixXd normal_check(2,1);
    normal_check << - std::sqrt(2.0) / 2.0,  - std::sqrt(2.0) / 2.0;
    EXPECT_TRUE(lf::dgfe::outwardNormal(a_edge).isApprox(normal_check));


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

TEST(integration, basic){

    // // (2^2) * (3^2) = 36
    // Eigen::MatrixXd a_point(2,1);
    // a_point <<      2, 3;
    // EXPECT_EQ(lf::dgfe::integrate(a_point, 2, 2), 36);

    
    //test simple triangle from paper
    Eigen::MatrixXd a_polygon(2,3);
    a_polygon <<    -1.0, 1.0, -1.0,
                    -1.0, 0.0, 1.0;
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 10), 0.0111339078, 1e-8);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 5), 0.0, 1e-8);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 20), -0.005890191, 1e-8);

    
}





} //namespace lf::dgfe::test