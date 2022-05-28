/**
 * @file
 * @brief Functionalities for reference basis functions of a discrete space in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/


#ifndef LEGENDRE_DGFE_H
#define LEGENDRE_DGFE_H

#include <lf/mesh/mesh.h>

#include <Eigen/Eigen>
#include <cmath>

#include "integration.h"
#include "bounding_box.h"



namespace lf::dgfe {

using scalar_t = double;
using size_type = lf::mesh::Mesh::size_type;

//container for legendre polynomial coefficients
//ith row indicates the ith polynomial
//jth column indicates the coefficient for x^j
const Eigen::Matrix<scalar_t, 3, 3> legendre_coeffs_{
    (Eigen::Matrix<scalar_t, 3, 3>() << 
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        -0.5,0.0,1.5).finished()};

/**
 * @brief returns the value of the ith legendre polynomial at x
 * 
 */
scalar_t legendre_polynomial(size_type i, scalar_t x);

/**
 * @brief returns the value of the definition (26) of the paper
 * https://link.springer.com/article/10.1007/s10915-018-0802-y
 * @param degree_p degree of the 1D polynomial space
 */
scalar_t C_i_j_k(size_type i, size_type j, size_type k, size_type degree_p);

/**
 * @brief returns a pair (degree_x, degree_y) derived from the basis function index of a cell
 * 
 * @param basis_index index of the basis fucntion on this cell
 * @param degree_of_space max degree of the space defined over the cell
 */
std::pair<size_type, size_type> multiIndexToDegree(size_type basis_index, size_type degree_of_space);

/**
 * @brief inverse function of multiIndexToDegree
 */
size_type degreeToMultiIndex(std::pair<size_type, size_type> degrees, size_type degree_of_space);


} //namespace lf::dgfe




#endif //LEGENDRE_DGFE_H