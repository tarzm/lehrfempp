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
#include "lf/mesh/test_utils/test_meshes.h"



namespace lf::dgfe::test {

using scalar_t = lf::dgfe::scalar_t;

TEST(legendre_dgfe, legendrePolynomials){
    scalar_t x = 0.8;
    scalar_t y = -0.7;
    EXPECT_EQ(legendre_polynomial(1, x) * legendre_polynomial(2,y), 0.18799999999999992);
}



} //namespace lf::dgfe:test