#ifndef LF_DGFEINTEGRATION_H
#define LF_DGFEINTEGRATION_H

/**
 * @file
 * @brief Functions for integrating polynomials on polygons
 * @author Tarzis Maurer
 * @date 2022-19-04
 * @copyright ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>

//degree of polynomials
using degree_t = unsigned;

namespace lf::dgfe {


template <typename SCALAR>
SCALAR integrate(Eigen::MatrixXd corners, degree_t degree_x, degree_t degree_y){
    auto n_cols = corners.cols();
    switch(n_cols){
        case 1: // entity is a point
            return std::pow(corners(0,0), degree_x) * std::pow(corners(1,0), degree_y);
            break;

        case 2: // entity is an edge
            SCALAR pre_factor = 1.0 / (1.0 + degree_x + degree_y);
            SCALAR sum = 0.0;
            sum += 
            for()
        

    }
}






} //namespace dgfe

#endif // LF_DGFEINTEGRATION_H