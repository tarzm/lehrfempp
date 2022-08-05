/**
 * @file
 * @brief Functionalities inlcuding average and jump operators
 * which are present in the general DGFEM form of linear degenerate
 * second order convection-diffusion-reaction bundary value problem
 * 
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef DGFE_AUX_OPERATORS_H
#define DGFE_AUX_OPERATORS_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>

#include "bounding_box.h"

namespace lf::dgfe {


/**
 * @brief computes the orthogonal L2 Projection of sqrt(a) * nabla_basis onto the discrete function space V^h
 * 
 * @tparam SCALAR scalar type of the function
 * @param elem_mat a placeholder for the local mass matrix
 * @param zeta_box_t the qr points used for integration mapped to the corresponding polygon
 * @param gram_dets_t gramian determinants used for quadrature rule
 * @param w_ref_t the weights of the quadrature rule
 * @param n_dofs number of local dofs
 * @param box bounding box of the polygon
 * @param a_evaluated the diffusion tensor evaluated at the qr points
 * @return LocalDofVectors two dof vectors which are the coefficients of the basis expansion in two dimensions
 */
template<typename SCALAR>
std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>> localL2ProjectionSqrtANablaBasis (Eigen::MatrixXd &zeta_box_t,
                                Eigen::VectorXd &gram_dets_t, Eigen::VectorXd &w_ref_t, unsigned &n_dofs,
                                lf::dgfe::BoundingBox &box, std::vector<Eigen::Matrix2d> a_evaluated){
    
    
    //setup of matrix and rhs
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>  elem_mat(n_dofs, n_dofs);
    elem_mat.setZero();
    //two rhs vectors for two dimensions
    Eigen::VectorXd rhs_0(elem_mat.cols());
    Eigen::VectorXd rhs_1(elem_mat.cols());
    rhs_0.setZero();
    rhs_1.setZero();

    unsigned max_legendre_degree = (n_dofs == 4) ? 1 : 2;

    //loop over trial basis functions
    for (int basis_trial = 0; basis_trial < n_dofs; basis_trial++){
        //loop over test basis functions
        for (int basis_test = 0; basis_test < n_dofs; basis_test++){
            //sum over qr points
            for (int i = 0; i < gram_dets_t.size(); i++){
                elem_mat(basis_trial, basis_test) +=    legendre_basis(basis_trial, max_legendre_degree, zeta_box_t.col(i)) 
                                                        * legendre_basis(basis_test, max_legendre_degree, zeta_box_t.col(i))
                                                        * w_ref_t[i] * gram_dets_t[i];
            }
        }
        //second loop over qr points for rhs assembly
        for(int i = 0; i < gram_dets_t.size(); i++){
            Eigen::Vector2d nabla_trial{legendre_basis_dx(basis_trial, max_legendre_degree, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                        legendre_basis_dy(basis_trial, max_legendre_degree, zeta_box_t.col(i)) * box.inverseJacobi(1)};

            rhs_0(basis_trial) += ((a_evaluated[i].sqrt()).row(0)).dot(nabla_trial) * w_ref_t[i] * gram_dets_t[i];
            rhs_1(basis_trial) += ((a_evaluated[i].sqrt()).row(1)).dot(nabla_trial) * w_ref_t[i] * gram_dets_t[i];
        }
    }

    Eigen::FullPivLU<decltype(elem_mat)> solver;
    solver.compute(elem_mat);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd phi_0 = solver.solve(rhs_0);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE 0 failed");
    Eigen::VectorXd phi_1 = solver.solve(rhs_1);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE 1 failed");

    return std::make_pair(phi_0, phi_1);
}

template<typename SCALAR>
std::vector<Eigen::Matrix<SCALAR, 2, 1>> evalL2Projection(std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>> &local_dof_vectors, Eigen::MatrixXd &local){
    auto dof_0 = local_dof_vectors.first;
    auto dof_1 = local_dof_vectors.second;

    unsigned max_legendre_degree = (dof_0.size() == 4) ? 1 : 2;

    std::vector<Eigen::Matrix<SCALAR, 2,1>> result(local.cols(), Eigen::Matrix<SCALAR, 2, 1>::Zero());

    //loop over basis functions
    for (int basis = 0; basis < dof_0.size(); basis++){
        //loop over local points
        for (int i = 0; i < local.cols(); i++){
            result[i][0] +=  dof_0[basis] * legendre_basis(basis, max_legendre_degree, local.col(i));
            result[i][1] +=  dof_1[basis] * legendre_basis(basis, max_legendre_degree, local.col(i));
        }
    }

    return result;
}

} //namespace lf::dgfe


#endif //DGFE_AUX_OPERATORS_H