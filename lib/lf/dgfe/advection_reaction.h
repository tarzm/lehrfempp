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

template<typename SCALAR, typename ADVECTION_COEFF, typename REACTION_COEFF, typename EDGESELECTOR>
class AdvectionReactionElementMatrixProvider {

public: 

    AdvectionReactionElementMatrixProvider(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, ADVECTION_COEFF b_coeff, REACTION_COEFF c_coeff,
                                            EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge,
                                            EDGESELECTOR boundary_minus_edge, unsigned integration_degree)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), b_coeff_(b_coeff), c_coeff_(c_coeff),
         boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         boundary_minus_edge_(std::move(boundary_minus_edge)) {
            LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
    }


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

        //local - global mapping
        lf::dgfe::BoundingBox box(cell);

        //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
        //     [nabla (b * w) + c*w] * v   over all cells

        //quadrule setup
        const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};

        //get sub-tessellation of cell
        auto sub_tessellation = subTessellation(&cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to global triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};
            //evaluation of the coefficient functions at qr points
            auto b = b_coeff_(cell, zeta_box_t);
            auto c = c_coeff_(cell, zeta_box_t);

            //loop over basis functions in trial space
            for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                //loop over basis functions in test space
                for(int basis_test = 0; basis_test < n_basis; basis_test++){
                    //sum over qr points
                    for (int i = 0; i < gram_dets_t.size(); i++){
                        //first part [nabla (b * w) + c*w] * v
                        elem_mat(basis_trial, basis_test) += ( b[i][0] * legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0)
                                                             + b[i][1] * legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)
                                                             + c[i] * legendre_basis(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) )
                                                             * legendre_basis(basis_test, max_legendre_degree_, zeta_box_t.col(i))
                                                             * w_ref_t[i] * gram_dets_t[i];
                    }
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

            // qr points mapped to gobal edge
            Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
            //gramian determinants
            Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};
            
            //coefficient evaluated at qr points
            auto b = b_coeff_(cell, zeta_box_s);

            //check wether weigthed sum of qr points satisfies
            //    (b(x) * n(x) < 0)
            //=> if so, edge belongs to the delta_minus_k set
            double sum = 0.0;
            for (int i = 0; i < gram_dets_s.size(); i++){
                sum += b[i].dot(normal) * w_ref_s[i] * gram_dets_s[i];
            }
            bool delta_minus_k = (sum < 0) ? true : false;

            // !!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!
            //    - ( b * n ) * upwind_jump[w] * v^+   over all edges which are not on boundary and belong to delta_minus_k

            if (!(boundary_edge_(*edge)) && delta_minus_k){    // edge must not be on boundary and belong to delta_minus_k

                //get pointer to polygon on other side of edge
                auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(edge);
                auto other_polygon = (polygon_pair.first.first == &cell) ? polygon_pair.second.first : polygon_pair.first.first;
                //get qr points mapped to other polygon's reference box
                lf::dgfe::BoundingBox box_other(*other_polygon);
                Eigen::MatrixXd zeta_box_other{box_other.inverseMap(zeta_global_s)};
                    
                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){

                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            elem_mat(basis_trial, basis_test) -=  (b[i].dot(normal))
                                                                * ( legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i))
                                                                    - legendre_basis(basis_trial, max_legendre_degree_, zeta_box_other.col(i)) )
                                                                * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                                                * w_ref_s[i] * gram_dets_s[i];
                        }
                    }
                }
            }
            // !!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!


            // !!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!
            //    - ( b * n ) * w^+ * v^+   over delta_minus_k and (boundary_d_edge or boundary_minus_edge)
            if (delta_minus_k && (boundary_d_edge_(*edge) || boundary_minus_edge_(*edge))){
                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){

                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            elem_mat(basis_trial, basis_test) -= (b[i].dot(normal))
                                                                * legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i))
                                                                * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                                                * w_ref_s[i] * gram_dets_s[i];
                    
                        }
                    }
                }
            }
            edge_sub_idx++;
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
    EDGESELECTOR boundary_edge_;
    EDGESELECTOR boundary_d_edge_;
    EDGESELECTOR boundary_minus_edge_;
};












} // namespace lf::dgfe

#endif // ADVECTION_REACTION_DGFE_H