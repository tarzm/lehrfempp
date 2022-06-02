/**
 * @file
 * @brief Tests of the the discontinuous galerkin finite elements
 * methods and their assembly alogrithms
 * @author Tarzis Maurer
 * @date May 22
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

using scalar_t = lf::dgfe::scalar_t;

TEST(legendre_dgfe, legendrePolynomials){
    scalar_t x = 0.8;
    scalar_t y = -0.7;
    EXPECT_EQ(legendre_polynomial(1, x) * legendre_polynomial(2,y), 0.18799999999999992);
}

TEST(legendre_dgfe, numberingBasisFunctions){
    //test implementation of getting the correct exponents in the degree O1 basis functions
    //from the multi index of the degrees
    std::vector<std::pair<int, int>> check_O1;
    check_O1.push_back(std::make_pair(0,0));
    check_O1.push_back(std::make_pair(0,1));
    check_O1.push_back(std::make_pair(1,0));
    check_O1.push_back(std::make_pair(1,1));

    for (int i = 0; i< 4; i++){
        auto degrees = multiIndexToDegree(i, 1);
        int basis_degree_x = degrees.first;
        int basis_degree_y = degrees.second;
        ASSERT_EQ(std::make_pair(basis_degree_x, basis_degree_y), check_O1[i]);
    }

    //test implementation of getting the correct exponents in the degree O2 basis functions
    //from the multi index of the degrees
    std::vector<std::pair<int, int>> check_O2;
    check_O2.push_back(std::make_pair(0,0));
    check_O2.push_back(std::make_pair(0,1));
    check_O2.push_back(std::make_pair(0,2));
    check_O2.push_back(std::make_pair(1,0));
    check_O2.push_back(std::make_pair(1,1));
    check_O2.push_back(std::make_pair(1,2));
    check_O2.push_back(std::make_pair(2,0));
    check_O2.push_back(std::make_pair(2,1));
    check_O2.push_back(std::make_pair(2,2));

    for (int i = 0; i< 9; i++){
        auto degrees = multiIndexToDegree(i, 2);
        int basis_degree_x = degrees.first;
        int basis_degree_y = degrees.second;
        ASSERT_EQ(std::make_pair(basis_degree_x, basis_degree_y), check_O2[i]);
    }

    //now check degreeToMultiIndex
    //O1
    for (int i = 0; i < 4; i++){
        auto index = degreeToMultiIndex(check_O1[i], 1);
        EXPECT_EQ(index, i);
    }
    //O2
    for (int i = 0; i < 9; i++){
        auto index = degreeToMultiIndex(check_O2[i], 2);
        EXPECT_EQ(index, i);
    }
}

// TEST(legendre_dgfe, C_i_j_k_O1){
//     //trial space basis
//     for (int i = 0; i < 4; i++){
//         //test space basis
//         for (int j = 0; j < 4; j++){
//             auto degrees_i = multiIndexToDegree(i, 1);
//             auto degrees_j = multiIndexToDegree(j, 1);
//             int i1 = degrees_i.first;
//             int i2 = degrees_i.second;
//             int j1 = degrees_j.first;
//             int j2 = degrees_j.second;
//             for (int k = 0; k <= i1 + j1; k++){
//                 for (int l = 0; l <= i2 + j2; l++){
//                     std::cout << "C_" << i1 << "_" << j1 << "_" << k << " = " <<  C_i_j_k(i1, j1, k, 1) << "\n";
//                 }
//             }
//         }
//     }
// }

} //namespace lf::dgfe:test