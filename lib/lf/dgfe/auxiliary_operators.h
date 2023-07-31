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

template<typename SCALAR, typename MESHFUNC>
struct L2ProjectionSqrtAGradBasis{
    std::vector<MESHFUNC> dim0;
    std::vector<MESHFUNC> dim1;

    std::vector<Eigen::Matrix<SCALAR, 2, 1>> operator()(const lf::mesh::Entity& e, const Eigen::MatrixXd& local, int basis) const {
        //loop over local points given
        std::vector<Eigen::Matrix<SCALAR, 2, 1>> result;
        auto vec_dim0 = dim0[basis](e, local);
        auto vec_dim1 = dim1[basis](e, local);
        for (int i = 0; i < local.cols(); i++){
            result.emplace_back(vec_dim0[i], vec_dim1[i]);
        }
        return result;
    }
};


template<typename SCALAR, typename MESH_FUNCTION, typename DGFE_MESHFUNCTION>
L2ProjectionSqrtAGradBasis<SCALAR, DGFE_MESHFUNCTION> l2_proj_diffusion(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, MESH_FUNCTION m_a_coeff, unsigned max_integration_degree){
    
    unsigned max_legendre_degree = dgfe_space_ptr->MaxLegendreDegree();
    unsigned n_local_dofs = (max_legendre_degree + 1) * (max_legendre_degree + 1);
    auto n_cells = dgfe_space_ptr->Mesh()->NumEntities(0);

    //DGFE dummy Meshfunction
    Eigen::VectorXd dummy_vec = Eigen::VectorXd::Zero(n_cells * n_local_dofs);
    lf::dgfe::MeshFunctionDGFE<SCALAR> dummy_MeshFunc(dgfe_space_ptr, dummy_vec);

    //setup of MeshFunctionDGFE vectors
    std::vector<decltype(dummy_MeshFunc)> dim_0_vec;
    std::vector<decltype(dummy_MeshFunc)> dim_1_vec;
    dim_0_vec.reserve(n_local_dofs);
    dim_1_vec.reserve(n_local_dofs);

    std::cout << "Originally dim_0_vec has length " << dim_0_vec.size() << "\n";

    //setup mass matrix
    unsigned n_dofs = dgfe_space_ptr->LocGlobMap().NumDofs();
    lf::assemble::COOMatrix<SCALAR> M(n_dofs, n_dofs);
    M.setZero();

    //quadrule setup
    const lf::quad::QuadRule qr_t = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), max_integration_degree);
    //qr points
    const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
    //weights
    Eigen::VectorXd w_ref_t{qr_t.Weights()};

    //dofhandler
    auto dofhandler = dgfe_space_ptr->LocGlobMap();
    
    //loop over cells
    for (const lf::mesh::Entity *cell : dgfe_space_ptr->Mesh()->Entities(0)){
        //dof info
        nonstd::span<const Eigen::Index> dof_plus = dofhandler.GlobalDofIndices(*cell);

        //local - global mapping
        lf::dgfe::BoundingBox box(*cell);
        //get sub-tessellation
        auto sub_tessellation = subTessellation(cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            //qr points mapped to triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

            //loop over basis functions in trial space
            for (int basis_trial = 0; basis_trial < n_local_dofs; basis_trial++){
                //loop over basis functions in test space
                for(int basis_test = 0; basis_test < n_local_dofs; basis_test++){
                    //sum over qr points
                    SCALAR sum = 0;
                    for (int i = 0; i < gram_dets_t.size(); i++){
                        sum += legendre_basis(basis_trial, max_legendre_degree, zeta_box_t.col(i))
                                * legendre_basis(basis_test, max_legendre_degree, zeta_box_t.col(i))
                                * w_ref_t[i] * gram_dets_t[i];
                    }
                    M.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum);
                }
            }
        }
    }

    ///////// MASS MATRIX IS ASSEMBLED
    //setup solver
    Eigen::SparseMatrix<SCALAR> M_crs = M.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<SCALAR>> solver;
    solver.compute(M_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");

    //rhs placeholder
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> rhs(n_dofs);

    //Loop over two dimensions
    for (int dim = 0; dim < 2; dim++){
        //loop over basis functions
        for (int basis = 0; basis < n_local_dofs; basis++){
            
            rhs.setZero();
            //loop over cells
            for (const lf::mesh::Entity *cell : dgfe_space_ptr->Mesh()->Entities(0)){
                //dof info
                nonstd::span<const Eigen::Index> dof_plus = dofhandler.GlobalDofIndices(*cell);

                //local - global mapping
                lf::dgfe::BoundingBox box(*cell);
                //get sub-tessellation
                auto sub_tessellation = subTessellation(cell);

                //loop over triangles in the sub-tessellation
                for(auto& tria_geo_ptr : sub_tessellation){
                    //qr points mapped to triangle
                    Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
                    // qr points mapped back into reference bounding box to retrieve values
                    Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
                    //gramian determinants
                    Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

                    //evaluated a
                    auto a_eval = m_a_coeff(*cell, zeta_box_t);

                    SCALAR sum = 0;
                    //sum over qr points
                    for (int i = 0; i < gram_dets_t.size(); i++){
                        Eigen::Vector2d nabla_basis{legendre_basis_dx(basis, max_legendre_degree, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                                    legendre_basis_dy(basis, max_legendre_degree, zeta_box_t.col(i)) * box.inverseJacobi(1)};

                        sum += (a_eval[i].sqrt() * nabla_basis)[dim] * w_ref_t[i] * gram_dets_t[i];
                    }
                    rhs[dof_plus[basis]] += sum;
                }
            }

            //solve LSE
            Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> sol_vec = solver.solve(rhs);
            LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

            //add mesh function to vector
            if (dim == 0){
                dim_0_vec.emplace_back(dgfe_space_ptr, sol_vec);
                std::cout << "Dim_0_vec has length " << dim_0_vec.size() << "\n";
            } else {
                dim_1_vec.emplace_back(dgfe_space_ptr, sol_vec);
            }
        }
    }

    L2ProjectionSqrtAGradBasis<SCALAR, decltype(dummy_MeshFunc)> result;
    result.dim0 = std::move(dim_0_vec);
    result.dim1 = std::move(dim_1_vec);

    std::cout << "Reslut L2 Projection has length dim0: " << result.dim0.size() << "\n";
    
    return result;
}




/**
 * @brief returns the L^2-projection of sqrt(a) * nabla v for all basis functions v
 */

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