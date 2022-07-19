/**
 * @file
 * @brief Diffusion element matrix provider
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef DIFFUSION_DGFE_H
#define DIFFUSION_DGFE_H

#include <lf/quad/quad.h>
#include <lf/mesh/utils/utils.h>

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"
#include "integration.h"
#include "mesh_function_dgfe.h"
#include "mesh_function_global.h"
#include "discontinuity_penalization.h"

namespace lf::dgfe {

template<typename SCALAR, typename DIFFUSION_COEFF, typename EDGESELECTOR>
class DiffusionElementMatrixProvider {

public: 

    DiffusionElementMatrixProvider(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, DIFFUSION_COEFF a_coeff,
                                    EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge, unsigned integration_degree,
                                    lf::dgfe::DiscontinuityPenalization disc_pen)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), a_coeff_(a_coeff),
         boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         disc_pen_(disc_pen), evaluated_edge_(std::move(initialize_evaluated_edge())) {
            LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
    }

    /**
     * @brief All cells are considered active in the default implementation
     */
    [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
        return true;
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> Eval(const lf::mesh::Entity &cell){

        const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
        //initialize element matrix
        Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(n_basis, n_basis);
        elem_mat.setZero();

        //local - global mapping
        lf::dgfe::BoundingBox box(cell);
        
        //quadrule setup
        const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};


        //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
        //     a * nabla(w) * nabla(v)  over all cells

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

            auto a = a_coeff_(cell, zeta_box_t);

            //loop over basis functions in trial space
            for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                //loop over bsis functions in test space
                for(int basis_test = 0; basis_test < n_basis; basis_test++){
                    //sum over qr points
                    for (int i = 0; i < gram_dets_t.size(); i++){

                        Eigen::Vector2d nabla_w{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                                legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};
                        Eigen::Vector2d nabla_v{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                                legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};

                        elem_mat(basis_trial, basis_test) += (a[i] * nabla_w).dot(nabla_v) * w_ref_t[i] * gram_dets_t[i];
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

        //loop over cells edges
        for (auto edge : cell.SubEntities(1)){
            if( (boundary_d_edge_(*edge) ||  !boundary_edge_(*edge)) && !evaluated_edge_(*edge) ){

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
            
                //get pointer to polygon on other side of edge
                auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(edge);
                auto other_polygon = (polygon_pair.first.first == &cell) ? polygon_pair.second.first : polygon_pair.first.first;
                //if edge is on boundary, other_polygon pointer is just to same cell
                if (other_polygon == nullptr){
                    other_polygon = &cell;
                }
                //get qr points mapped to other polygon's reference box
                lf::dgfe::BoundingBox box_other(*other_polygon);
                Eigen::MatrixXd zeta_box_other{box_other.inverseMap(zeta_global_s)};
                Eigen::Vector2d normal_other = normal *= -1.0;

                //calculate A_F
                Eigen::MatrixXd A_F_mat = Eigen::MatrixXd::Zero(2, gram_dets_s.size());
                for (int i = 0; i < gram_dets_s.size(); i++){
                    A_F_mat.col(i) = a_coeff_(cell, zeta_box_s.col(i))[0].sqrt() * normal;
                }
                SCALAR A_F = A_F_mat.lpNorm<Eigen::Infinity>();
                A_F = A_F * A_F;

                
                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){

                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){

                            Eigen::Vector2d jump_trial{legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * normal
                                                     + legendre_basis(basis_trial, max_legendre_degree_, zeta_box_other.col(i)) * normal_other};
                            
                            Eigen::Vector2d jump_test{legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * normal
                                                     + legendre_basis(basis_test, max_legendre_degree_, zeta_box_other.col(i)) * normal_other};

                            //if edge is on edge_d, the jump operators need to be divided by 2
                            //because jump operators are evaluated just like there was a polygon on other side of edge
                            if(other_polygon == &cell){
                                jump_trial *= 0.5;
                                jump_test *= 0.5;
                            }

                            //!!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!!
                            elem_mat(basis_trial, basis_test) += disc_pen_(*edge, A_F) * jump_trial.dot(jump_test) * w_ref_s[i] * gram_dets_s[i];
                            //!!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!!
                            
                            //!!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!!
                            Eigen::Vector2d nabla_trial{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                        legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            Eigen::Vector2d nabla_test{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                       legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            Eigen::Vector2d nabla_trial_other{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_other.col(i)) * box_other.inverseJacobi(0),
                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_other.col(i)) * box_other.inverseJacobi(1)};

                            Eigen::Vector2d nabla_test_other{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_other.col(i)) * box_other.inverseJacobi(0),
                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_other.col(i)) * box_other.inverseJacobi(1)};
                                                            
                            Eigen::Vector2d average_a_nabla_trial{0.5 * (a_coeff_(cell, zeta_box_s.col(i))[0] * (nabla_trial + nabla_trial_other))};

                            Eigen::Vector2d average_a_nabla_test{0.5 * (a_coeff_(cell, zeta_box_s.col(i))[0] * (nabla_test + nabla_test_other))};

                            elem_mat(basis_trial, basis_test) -= (average_a_nabla_trial.dot(jump_test) * average_a_nabla_test.dot(jump_trial)) * w_ref_s[i] * gram_dets_s[i];
                            //!!!!!!!!!!!!! END THIRD TERM !!!!!!!!!!!!!!
                        }
                    }
                }
                evaluated_edge_(*edge) = true;
            } // second and third term if
            edge_sub_idx++;
        }
        return elem_mat;
    }

private:

    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
    unsigned integration_degree_;
    unsigned max_legendre_degree_;
    lf::quad::QuadRuleCache qr_cache_;
    DIFFUSION_COEFF a_coeff_;
    EDGESELECTOR boundary_edge_;
    EDGESELECTOR boundary_d_edge_;
    lf::dgfe::DiscontinuityPenalization disc_pen_;
    //used to make sure the second and third term is only evaluated once per edge
    //is true if edge has already been evaluated
    EDGESELECTOR evaluated_edge_;

    lf::mesh::utils::CodimMeshDataSet<bool> initialize_evaluated_edge(){
        lf::mesh::utils::CodimMeshDataSet<bool> result(dgfe_space_ptr_->Mesh(), 1, false);
        return result;
    }

};

} //namespace lf::dgfe

#endif // DIFFUSION_DGFE_H