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
#include <lf/assemble/assemble.h>

#include "bounding_box.h"
#include "dgfe_providers.h"
#include "mesh_function_dgfe.h"

namespace lf::dgfe {

template<typename SCALAR, typename MESH_FUNCTION>
std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>, std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>> L2ProjectionSqrtANablaBasis(
            std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, MESH_FUNCTION m_a_coeff, unsigned max_integration_degree){
    
    unsigned max_legendre_degree = dgfe_space_ptr->MaxLegendreDegree();
    unsigned n_local_dofs = (max_legendre_degree + 1) * (max_legendre_degree + 1);

    //setup of MeshFunctionDGFE vectors
    std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>> dim_0_vec;
    std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>> dim_1_vec;

    //setup mass matrix
    unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
    lf::assemble::COOMatrix<SCALAR> M(n_dofs, n_dofs);
    M.setZero();
    lf::dgfe::DGFEMassElementMatrixST<SCALAR> massMatrixProvider(max_integration_degree, max_legendre_degree);
    lf::assemble::AssembleMatrixLocally(0, dgfe_space_ptr->LocGlobMap(), dgfe_space_ptr->LocGlobMap(), massMatrixProvider, M);

    //setup solver
    Eigen::SparseMatrix<SCALAR> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<SCALAR>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");

    //rhs placeholder
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> rhs(n_dofs);
    
    //loop over dimensions
    for (int dim = 0; dim < 2; dim++){
        //loop over basis functions
        for (int basis = 0; basis < n_local_dofs; basis++){

            rhs.setZero();
            //assemble rhs
            lf::dgfe::L2ProjectionSqrtANablaBasisLoadVector<SCALAR, decltype(m_a_coeff)> l2_projection_provider(dgfe_space_ptr, m_a_coeff, dim, basis, max_integration_degree);
            lf::assemble::AssembleVectorLocally(0, dgfe_space_ptr->LocGlobMap(), l2_projection_provider, rhs);
            //solve LSE
            Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> sol_vec = solver.solve(rhs);
            LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

            //fill MeshFunction vector
            if(dim == 0){
                dim_0_vec.emplace_back(dgfe_space_ptr, sol_vec);
            } else {
                dim_1_vec.emplace_back(dgfe_space_ptr, sol_vec);
            }
        }
    }

    return std::make_pair(dim_0_vec, dim_1_vec);
}

} //namespace lf::dgfe


#endif //DGFE_AUX_OPERATORS_H