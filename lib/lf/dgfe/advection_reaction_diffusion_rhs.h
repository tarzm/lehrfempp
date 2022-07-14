/**
 * @file
 * @brief Advection reaction diffusion element vector provider
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/


#ifndef ADVECTION_REACTION_DIFFUSION_RHS_H
#define ADVECTION_REACTION_DIFFUSION_RHS_H

#include <lf/quad/quad.h>
#include <lf/mesh/utils/utils.h>

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"
#include "integration.h"
#include "mesh_function_dgfe.h"
#include "mesh_function_global.h"

namespace lf::dgfe {

template<typename SCALAR, typename DIFFUSION_COEFF, typename ADVECTION_COEFF, typename EDGESELECTOR, typename FUNCTOR_F, typename FUNCTOR_G_D, typename FUNCTOR_G_N>
class AdvectionReactionDiffusionRHS {

public:
    AdvectionReactionDiffusionRHS(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, FUNCTOR f, FUNCTOR_G_D gD, FUNCTOR_G_N gN,
                                    DIFFUSION_COEFF a_coeff, ADVECTION_COEFF b_coeff, 
                                    EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge,
                                    EDGESELECTOR boundary_n_edge, unsigned integration_degree)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), b_coeff_(b_coeff), a_coeff_(a_coeff),
         boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         boundary_n_edge_(std::move(boundary_n_edge)), f_(f), gD_(gD), gN_(gN) {
            LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
    }


    /**
     * @brief All cells are considered active in the default implementation
     */
    [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
        return true;
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> Eval(const lf::mesh::Entity &cell) const {
        const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
        //initialize element vector
        Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> elem_vec(n_basis);
        elem_mat.setZero();

        //local - global mapping
        lf::dgfe::BoundingBox box(cell);


        //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
        //     f * w   over all cells

        //quadrule setup
        const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

            auto f_evaluated = f_(cell, zeta_box_t);

            //loop over basis functions in test space
            for (int basis_test = 0; basis_test < n_basis; basis_test++){
                
                //sum over qr points
                for (int i = 0; i < gram_dets_t.size(); i++){
                    //first part [nabla (b * w) + c*w] * v
                    elem_mat(basis_test) += * w_ref_t[i] * gram_dets_t[i] * f_evaluated[i] * legendre_basis(basis_test, max_legendre_degree_, zeta_box_t.col(i));
                }
            }
        }
        //!!!!!!!!!!!!! END FIRST TERM !!!!!!!!!!!!!!

        //quadrule setup
        const lf::quad::QuadRule qr_s = qr_cache_.Get(lf::base::RefEl::kSegment(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_s{qr_s.Points()};
        //weights
        Eigen::VectorXd w_ref_s{qr_s.Weights()};

        size_type edge_sub_idx = 0;

        //loop over cell's edges
        for (auto edge : cell.SubEntities(1)){

            //normal n
            auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
            //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
            normal *= (int) (cell.RelativeOrientations()[edge_sub_idx]);

            // qr points mapped to segment
            Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
            //gramian determinants
            Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};

            auto b = b_coeff_(cell, zeta_box_s);
            auto gD_evaluated = gD_(cell, zeta_box_s);

            //check wether weigthed sum of qr points satisfies
            //    (b(x) * n(x) < 0)
            //=> if so, edge belongs to the delta_minus_k set of that polygon
            double sum;
            for (int i = 0; i < gram_dets_s.size(); i++){
                sum += b[i].dot(normal) * w_ref_s[i] * gram_dets_s[i];
            }
            bool delta_minus_k = (sum < 0) ? true : false;

            // !!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!
            //    - ( b * n ) * g_D * v^+   over all edges which are either on boundary_D or boundary_minus and belong to delta_minus_k

            if (delta_minus_k && (boundary_d_edge(*edge) || (boundary_minus_edge(*edge)))){

                //loop over bsis functions in test space
                for(int basis_test = 0; basis_test < n_basis; basis_test++){
                    //sum over qr points
                    for (int i = 0; i < gram_dets_s.size(); i++){
                        elem_vec -= (b[i][0] * normal[0] + b[i][1] * normal[1]) * gD_evaluated[i]
                                     * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                     * w_ref_s[i] * gram_dets_s[i];
                    }
                }
            }
            edge_sub_idx++;
        }
        return elem_vec;
    } 


private:
    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
    unsigned integration_degree_;
    unsigned max_legendre_degree_;
    lf::quad::QuadRuleCache qr_cache_;
    FUNCTOR f_;
    FUNCTOR_G_D gD_;
    FUNCTOR_G_N gN_;
    DIFFUSION_COEFF a_coeff_;
    ADVECTION_COEFF b_coeff_;
    EDGESELECTOR boundary_edge_;
    EDGESELECTOR boundary_d_edge_;
    EDGESELECTOR boundary_n_edge_;
};


















} //namespace lf::dgfe

#define ADVECTION_REACTION_DIFFUSION_RHS_H