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
#include "integration.h"

namespace lf::dgfe {

//polyomial expansion of a function in 2D
//format: {[coefficient, (degree x, degree y)], [coefficient, (degree x, degree y)], ... }
using Polynomial = std::vector<std::pair<scalar_t, std::pair<size_type, size_type>>>;


/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/dgfe/dgfe_providers.h
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

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/dgfe/dgfe_providers.h
 * @brief Computing the element matrix for the mass matrix using the SubTessellationIntegrator class
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref AssembleMatrixLocally().
 */
class DGFEMassElementMatrixST {
    public:
        /**
         * @brief Construct a new DGFEMassElementMatrixST object
         * 
         * @param integration_max_degree the maximum degree used for the quadrature rule
         * @param legendre_max_degree maximum degree of legendre polynomials: Either 1 or 2
         */
        DGFEMassElementMatrixST(unsigned max_integration_degree, unsigned max_legendre_degree);

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
         * @return a 4x4 or a 9x9 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> Eval(const lf::mesh::Entity &cell) const;
    
    private:

        unsigned max_integration_degree_;
        unsigned max_legendre_degree_;
};

} //namespace lf::dgfe


#endif //DGFE_PROVIDERS_H