/**
 * @file
 * @brief Advection reaction element matrix provider
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef ADVECTION_REACTION_DGFE_H
#define ADVECTION_REACTION_DGFE_H

#include <lf/quad/quad.h>
#include <lf/mesh/utils/utils.h>

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"
#include "integration.h"
#include "mesh_function_dgfe.h"
#include "mesh_function_global.h"

namespace lf::dgfe {

template<typename SCALAR, typename ADVECTION_COEFF, typename REACTION_COEFF>
class AdvectionReactionElementMatrixProvider {

public: 

    AdvectionReactionElementMatrixProvider(std::shared_ptr<const lf::dgfe::DGFESPace> dgfe_space_ptr, ADVECTION_COEFF b_coeff, REACTION_COEFF c_coeff, unsigned integration_degree)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree), max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()) {}


    /**
     * @brief All cells are considered active in the default implementation
     */
    [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
        return true;
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> Eval(const lf::mesh::Entity &cell) const{

        unsigned matrix_size = (max_legendre_degree_ == 1) ? 4 : 9;
        Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(matrix_size, matrix_size);

        int i1;
        int i2;
        int j1;
        int j2;
        auto eval_lambda = [&i1, &i2, &j1, &j2](const lf::mesh::Entity &entity, Eigen::MatrixXd &local) -> std::vector<SCALAR> {
            std::vector<SCALAR> result(local.cols());
            for (int i = 0; i < local.cols(); i++){
                result.at(i) = lf::dgfe::legendre_polynomial_2D(i1, i2, local.col(i)) * lf::dgfe::legendre_polynomial_2D(j1, j2, local.col(i));
            }
            return result;
        };


    }




private:
    std::shared_ptr<const lf::dgfe::DGFESPace> dgfe_space_ptr_;
    unsigned integration_degree_;
    unsigned max_legendre_degree_;


};












}

#endif // ADVECTION_REACTION_DGFE_H