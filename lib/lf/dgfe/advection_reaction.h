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

    AdvectionReactionElementMatrixProvider(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, ADVECTION_COEFF b_coeff, REACTION_COEFF c_coeff, unsigned integration_degree)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree), max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), b_coeff_(b_coeff), c_coeff_(c_coeff) {}


    /**
     * @brief All cells are considered active in the default implementation
     */
    [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
        return true;
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> Eval(const lf::mesh::Entity &cell) const{

        const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
        //initialize element matrix
        Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(n_basis, n_basis);
        elem_mat.setZero();

        //quadrule setup
        const lf::quad::QuadRule qr = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);

        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);
        // qr points
        const Eigen::MatrixXd zeta_ref{qr.Points()};
        //weights
        Eigen::VectorXd w_ref{qr.Weights()};
        //local - global mapping
        lf::dgfe::BoundingBox box(cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global{tria_geo_ptr->Global(zeta_ref)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box{box.inverseMap(zeta_global)};
            //gramian determinants
            Eigen::VectorXd gram_dets{tria_geo_ptr->IntegrationElement(zeta_ref)};

            auto b = b_coeff_(cell, zeta_box);
            auto c = c_coeff_(cell, zeta_box);
            //loop over basis functions in trial space
            for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                //loop over bsis functions in test space
                for(int basis_test = 0; basis_test < n_basis; basis_test++){

                    //sum over qr points
                    for (int i = 0; i < qr.Points().cols(); i++){
                        //first part [nabla (b * w) + c*w] * v
                        elem_mat(basis_trial, basis_test) += ( b[i][0] * legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box.col(i))
                                                             + b[i][1] * legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box.col(i))
                                                             + c[i] * legendre_basis(basis_trial, max_legendre_degree_, zeta_box.col(i)) )
                                                             * legendre_basis(basis_test, max_legendre_degree_, zeta_box.col(i));
                    }

                }
            }
        }
        return elem_mat;
    }




private:
    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
    unsigned integration_degree_;
    unsigned max_legendre_degree_;
    lf::quad::QuadRuleCache qr_cache_;
    ADVECTION_COEFF b_coeff_;
    REACTION_COEFF c_coeff_;

};












}

#endif // ADVECTION_REACTION_DGFE_H