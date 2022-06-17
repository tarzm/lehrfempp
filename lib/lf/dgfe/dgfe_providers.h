/**
 * @file
 * @brief Classes supporting assembly of element matrices/vectors in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#ifndef DGFE_PROVIDERS_H
#define DGFE_PROVIDERS_H

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"

namespace lf::dgfe {

//polyomial expansion of a function in 2D
//format: {[coefficient, (degree x, degree y)], [coefficient, (degree x, degree y)], ... }
using Polynomial = std::vector<std::pair<scalar_t, std::pair<size_type, size_type>>>;


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
        [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        /**
         * @brief main routine for the computation of element matrices
         *
         * @param cell reference to the polytopic cell for
         *        which the element matrix should be computed.
         * @return a 4x4 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, 4, 4> Eval(const lf::mesh::Entity &cell) const;
};

class DGFEO1LocalLoadVector {
    public:
        using ElemVec = Eigen::Matrix<scalar_t, 4, 1>;

        DGFEO1LocalLoadVector(Polynomial polynomial) : polynomial_(polynomial) {}

        bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        ElemVec Eval(const lf::mesh::Entity &cell) const;

    private:
        Polynomial polynomial_;


};

class DGFEO2MassElementMatrix {
    public:
        /**
         * @brief All cells are considered active in the default implementation
         */
        [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        /**
         * @brief main routine for the computation of element matrices
         *
         * @param cell reference to the polytopic cell for
         *        which the element matrix should be computed.
         * @return a 9x9 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, 9, 9> Eval(const lf::mesh::Entity &cell) const;
};

class DGFEO2LocalLoadVector {
    public:
        using ElemVec = Eigen::Matrix<scalar_t, 9, 1>;
        

        DGFEO2LocalLoadVector(Polynomial polynomial) : polynomial_(polynomial) {}

        bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        ElemVec Eval(const lf::mesh::Entity &cell) const;

    private:
        Polynomial polynomial_;


};


} //namespace lf::dgfe


#endif //DGFE_PROVIDERS_H