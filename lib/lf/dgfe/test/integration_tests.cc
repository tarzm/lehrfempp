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

TEST(integration, basic){

    // (2^2) * (3^2) = 36
    Eigen::MatrixXd a_point(2,1);
    a_point <<      2, 3;
    EXPECT_EQ(lf::dgfe::integrate(a_point, 2, 2), 36);

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




    
}





} //namespace lf::dgfe::test