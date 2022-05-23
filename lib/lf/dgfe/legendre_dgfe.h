/**
 * @file
 * @brief Classes supporting the discontinuous finite element discretization of boundary value problems
 * with mapped legendre basis functions
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/


#ifndef LEGENDRE_DGFE_H
#define LEGENDRE_DGFE_H

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <cmath>

#include "integration.h"
#include "bounding_box.h"

#include <lf/mesh/mesh.h>

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
 * @ingroup entity_matrix_provider
 * @headerfile lf/dgfe/legendre_dgfe.h
 * @brief Computing the element matrix for the mass matrix
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref AssembleMatrixLocally().
 */
class DGFEO1MassElementMatrix {
    public:
        /**
         * @brief All cells are considered active in the default implementation
         */
        [[nodiscard]] virtual bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        /**
         * @brief main routine for the computation of element matrices
         *
         * @param cell reference to the polytopic cell for
         *        which the element matrix should be computed.
         * @return a 2x2 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, 4, 4> Eval(const lf::mesh::Entity &cell) const;
};



} //namespace lf::dgfe




#endif //LEGENDRE_DGFE_H