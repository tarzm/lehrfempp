
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
#include "auxiliary_operators.h"

namespace lf::dgfe {

template<typename SCALAR, typename DIFFUSION_COEFF, typename EDGESELECTOR>
class DiffusionElementMatrixProvider {

public: 

    using l2_proj_sqrt_a_nabla_basis = std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>, std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>>;

    DiffusionElementMatrixProvider(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, DIFFUSION_COEFF a_coeff,
                                    EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge, unsigned integration_degree,
                                    lf::dgfe::DiscontinuityPenalization disc_pen, l2_proj_sqrt_a_nabla_basis &l2_proj)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), a_coeff_(a_coeff),
         boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         disc_pen_(disc_pen), evaluated_edge_(std::move(initialize_evaluated_edge())), l2_projection_(l2_proj) {
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
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};
        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);


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

            //diffusion coefficient matrix evaluated at qr points
            auto a = a_coeff_(cell, zeta_box_t);

            //loop over basis functions in trial space
            for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                //loop over basis functions in test space
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
    l2_proj_sqrt_a_nabla_basis &l2_projection_;


    lf::mesh::utils::CodimMeshDataSet<bool> initialize_evaluated_edge(){
        lf::mesh::utils::CodimMeshDataSet<bool> result(dgfe_space_ptr_->Mesh(), 1, false);
        return result;
    }

};

template<typename TMPMATRIX, typename SCALAR, typename DIFFUSION_COEFF, typename EDGESELECTOR>
class DiffusionMatrixAssembler{

public:
    using l2_proj_sqrt_a_nabla_basis = std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>, std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>>;

    DiffusionMatrixAssembler(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, DIFFUSION_COEFF a_coeff,
                                    EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge, unsigned integration_degree,
                                    lf::dgfe::DiscontinuityPenalization disc_pen, l2_proj_sqrt_a_nabla_basis l2_proj)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), a_coeff_(a_coeff),
         boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         disc_pen_(disc_pen), evaluated_edge_(std::move(initialize_evaluated_edge())), l2_projection_(l2_proj) {
            LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
    }

    void assemble(TMPMATRIX &matrix){
        //FIRST PART => LOOP OVER CELLS FOR VOLUME PART
        //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
        //     a * nabla(w) * nabla(v)  over all cells
        const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
        //initialize element matrix
        Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(n_basis, n_basis);

        //quadrule setup
        const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};

        auto dofhandler = dgfe_space_ptr_->LocGlobMap();
        
        //loop over cells
        for (const lf::mesh::Entity *cell : dgfe_space_ptr_->Mesh()->Entities(0)){
            
            //reset element matrix
            elem_mat.setZero();

            //local - global mapping
            lf::dgfe::BoundingBox box(*cell);
            //get sub-tessellation
            auto sub_tessellation = subTessellation(cell);
            
            //loop over triangles in the sub-tessellation
            for(auto& tria_geo_ptr : sub_tessellation){
                // qr points mapped to triangle
                Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
                //gramian determinants
                Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

                //diffusion coefficient matrix evaluated at qr points
                auto a = a_coeff_(*cell, zeta_box_t);

                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over basis functions in test space
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
            // row indices of for contributions of cells
            nonstd::span<const Eigen::Index> row_idx(dofhandler.GlobalDofIndices(*cell));
            // Column indices of for contributions of cells
            nonstd::span<const Eigen::Index> col_idx(dofhandler.GlobalDofIndices(*cell));
            //assembly double loop
            for (int i = 0; i < n_basis; i++) {
                for (int j = 0; j < n_basis; j++) {
                // Add the element at position (i,j) of the local matrix
                // to the entry at (row_idx[i], col_idx[j]) of the global matrix
                matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));
                }
            }  // end assembly local double loop
        }
        //!!!!!!!!!!!!! END FIRST TERM !!!!!!!!!!!!!!

        //quadrule setup
        const lf::quad::QuadRule qr_s = qr_cache_.Get(lf::base::RefEl::kSegment(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_s{qr_s.Points()};
        //weights
        Eigen::VectorXd w_ref_s{qr_s.Weights()};

        //loop over edges
        for (const lf::mesh::Entity *edge : dgfe_space_ptr_->Mesh()->Entities(1)){

            auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(edge);
            const lf::mesh::Entity &cell = *polygon_pair.first.first;

            //normal n
            auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
            //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
            normal *= (int) (cell.RelativeOrientations()[polygon_pair.first.second]);

            lf::dgfe::BoundingBox box(cell);
            // qr points mapped to segment
            Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
            //gramian determinants
            Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};

            //calculate A_F
            Eigen::MatrixXd A_F_mat = Eigen::MatrixXd::Zero(2, gram_dets_s.size());
            for (int i = 0; i < gram_dets_s.size(); i++){
                A_F_mat.col(i) = a_coeff_(cell, zeta_box_s.col(i))[0].sqrt() * normal;
            }
            SCALAR A_F = A_F_mat.lpNorm<Eigen::Infinity>();
            A_F = A_F * A_F;

            //discontinuity penalization
            auto disc_pen = disc_pen_(*edge, A_F);


            //!!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!!
            if (!boundary_edge_(*edge)){ //interior edge
                
                //pointers to adjacent polygons and the edge's local index in those
                auto polygon_plus = polygon_pair.first.first;
                auto polygon_minus = polygon_pair.second.first;
                auto edge_sub_idx_plus = polygon_pair.first.second;
                auto edge_sub_idx_minus = polygon_pair.second.second;

                //dof info
                nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                nonstd::span<const Eigen::Index> dof_minus(dofhandler.GlobalDofIndices(*polygon_minus));

                //local - global mappings
                lf::dgfe::BoundingBox box_plus(*polygon_plus);
                lf::dgfe::BoundingBox box_minus(*polygon_minus);

                //mapped qr points to respective bounding box
                Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
                Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};

                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over basis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){
                        
                        SCALAR sum = 0.0;
                        //First part wi+ * vi+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum * disc_pen);
                        
                        sum = 0.0;
                        //Second part - wi+ * vj+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -sum * disc_pen);

                        sum = 0.0;
                        //Third part - wj+ * vi+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -sum * disc_pen);

                        sum = 0.0;
                        //Third part + wj+ * vj+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], sum * disc_pen);
                    }
                }
            }
            //still second term but for boundary edges => simpler expression
            if (boundary_d_edge_(*edge)){

                auto polygon_plus = polygon_pair.first.first;
                //dof info
                nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                
                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){

                        SCALAR sum = 0.0;
                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum * disc_pen);
                    }
                }
            }
            //!!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!!


            //!!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!!

            //diffusion coefficient matrix evaluated at qr points
            auto a = a_coeff_(cell, zeta_box_s);

            if (!boundary_edge_(*edge)){
                
                //pointers to adjacent polygons and the edge's local index in those
                auto polygon_plus = polygon_pair.first.first;
                auto polygon_minus = polygon_pair.second.first;
                auto edge_sub_idx_plus = polygon_pair.first.second;
                auto edge_sub_idx_minus = polygon_pair.second.second;

                //dof info
                nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                nonstd::span<const Eigen::Index> dof_minus(dofhandler.GlobalDofIndices(*polygon_minus));

                //local - global mappings
                lf::dgfe::BoundingBox box_plus(*polygon_plus);
                lf::dgfe::BoundingBox box_minus(*polygon_minus);

                //mapped qr points to respective bounding box
                Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
                Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};

                auto normal_plus = normal;
                auto normal_minus = -1.0 * normal_plus;

                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){
                        
                        //First, test functions on plus side are nonzero
                        ///////////////////////////////////////////////
                        SCALAR sum = 0.0;
    
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -0.5 * sum);

                        //Secondly, test functions on minus side are nonzero
                        ////////////////////////////////////////////////////
                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], -0.5 * sum);
                    }
                }
            }

            //still third term but for boundary edges => simpler expression
            if (boundary_d_edge_(*edge)){

                auto polygon_plus = polygon_pair.first.first;
                //dof info
                nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                
                //loop over basis functions in trial space
                for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(int basis_test = 0; basis_test < n_basis; basis_test++){

                        SCALAR sum = 0.0;
                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            sum +=      (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * normal)
                                    +   (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * normal)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -sum);
                    }
                }
            }
                


        }


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
    l2_proj_sqrt_a_nabla_basis &l2_projection_;

    lf::mesh::utils::CodimMeshDataSet<bool> initialize_evaluated_edge(){
        lf::mesh::utils::CodimMeshDataSet<bool> result(dgfe_space_ptr_->Mesh(), 1, false);
        return result;
    }
};

} //namespace lf::dgfe

#endif // DIFFUSION_DGFE_H 


// /**
//  * @file
//  * @brief Diffusion element matrix provider
//  * @author Tarzis Maurer
//  * @date June 22
//  * @copyright ETH Zurich
// */

// #ifndef DIFFUSION_DGFE_H
// #define DIFFUSION_DGFE_H

// #include <lf/quad/quad.h>
// #include <lf/mesh/utils/utils.h>

// #include "legendre_dgfe.h"
// #include "dgfe_space.h"
// #include "bounding_box.h"
// #include "integration.h"
// #include "mesh_function_dgfe.h"
// #include "mesh_function_global.h"
// #include "discontinuity_penalization.h"
// #include "auxiliary_operators.h"

// namespace lf::dgfe {

// template<typename TMPMATRIX, typename SCALAR, typename DIFFUSION_COEFF, typename EDGESELECTOR>
// class DiffusionMatrixAssembler{

// public:
//     using l2_proj_sqrt_a_nabla_basis = std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>, std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>>;

//     DiffusionMatrixAssembler(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, DIFFUSION_COEFF a_coeff,
//                                     EDGESELECTOR boundary_edge, EDGESELECTOR boundary_d_edge, unsigned integration_degree,
//                                     lf::dgfe::DiscontinuityPenalization disc_pen, l2_proj_sqrt_a_nabla_basis l2_proj)
//         : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
//          max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), a_coeff_(a_coeff),
//          boundary_edge_(std::move(boundary_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
//          disc_pen_(disc_pen), evaluated_edge_(std::move(initialize_evaluated_edge())), l2_projection_(l2_proj) {
//             LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
//     }

//     void assemble(TMPMATRIX &matrix){
//         //FIRST PART => LOOP OVER CELLS FOR VOLUME PART
//         //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
//         //     a * nabla(w) * nabla(v)  over all cells
//         const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
//         //initialize element matrix
//         Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(n_basis, n_basis);

//         //quadrule setup
//         const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
//         // qr points
//         const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
//         //weights
//         Eigen::VectorXd w_ref_t{qr_t.Weights()};

//         auto dofhandler = dgfe_space_ptr_->LocGlobMap();
        
//         //loop over cells
//         for (const lf::mesh::Entity *cell : dgfe_space_ptr_->Mesh()->Entities(0)){
            
//             //reset element matrix
//             elem_mat.setZero();

//             //local - global mapping
//             lf::dgfe::BoundingBox box(*cell);
//             //get sub-tessellation
//             auto sub_tessellation = subTessellation(cell);
            
//             //loop over triangles in the sub-tessellation
//             for(auto& tria_geo_ptr : sub_tessellation){
//                 // qr points mapped to triangle
//                 Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
//                 // qr points mapped back into reference bounding box to retrieve values
//                 Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
//                 //gramian determinants
//                 Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

//                 //diffusion coefficient matrix evaluated at qr points
//                 auto a = a_coeff_(*cell, zeta_box_t);

//                 //loop over basis functions in trial space
//                 for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
//                     //loop over basis functions in test space
//                     for(int basis_test = 0; basis_test < n_basis; basis_test++){
//                         //sum over qr points
//                         for (int i = 0; i < gram_dets_t.size(); i++){

//                             Eigen::Vector2d nabla_w{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
//                                                     legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};
//                             Eigen::Vector2d nabla_v{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
//                                                     legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};

//                             elem_mat(basis_trial, basis_test) += (a[i] * nabla_w).dot(nabla_v) * w_ref_t[i] * gram_dets_t[i];
//                         }
//                     }
//                 }
//             }
//             // row indices of for contributions of cells
//             nonstd::span<const Eigen::Index> row_idx(dofhandler.GlobalDofIndices(*cell));
//             // Column indices of for contributions of cells
//             nonstd::span<const Eigen::Index> col_idx(dofhandler.GlobalDofIndices(*cell));
//             //assembly double loop
//             for (int i = 0; i < n_basis; i++) {
//                 for (int j = 0; j < n_basis; j++) {
//                 // Add the element at position (i,j) of the local matrix
//                 // to the entry at (row_idx[i], col_idx[j]) of the global matrix
//                 matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));
//                 }
//             }  // end assembly local double loop
//         }
//         //!!!!!!!!!!!!! END FIRST TERM !!!!!!!!!!!!!!

//         //quadrule setup
//         const lf::quad::QuadRule qr_s = qr_cache_.Get(lf::base::RefEl::kSegment(), integration_degree_);
//         // qr points
//         const Eigen::MatrixXd zeta_ref_s{qr_s.Points()};
//         //weights
//         Eigen::VectorXd w_ref_s{qr_s.Weights()};

//         //loop over edges
//         for (const lf::mesh::Entity *edge : dgfe_space_ptr_->Mesh()->Entities(1)){

//             auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(edge);
//             const lf::mesh::Entity &cell = *polygon_pair.first.first;

//             //normal n
//             auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
//             //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
//             normal *= (int) (cell.RelativeOrientations()[polygon_pair.first.second]);

//             lf::dgfe::BoundingBox box(cell);
//             // qr points mapped to segment
//             Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
//             // qr points mapped back into reference bounding box to retrieve values
//             Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
//             //gramian determinants
//             Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};

//             //calculate A_F
//             Eigen::MatrixXd A_F_mat = Eigen::MatrixXd::Zero(2, gram_dets_s.size());
//             for (int i = 0; i < gram_dets_s.size(); i++){
//                 A_F_mat.col(i) = a_coeff_(cell, zeta_box_s.col(i))[0].sqrt() * normal;
//             }
//             SCALAR A_F = A_F_mat.lpNorm<Eigen::Infinity>();
//             A_F = A_F * A_F;

//             ///////calculate AF new
//             //get midpoint of edge
//             // auto corners = lf::geometry::Corners(*(edge->Geometry()));
//             // auto midpoint = (corners.col(0) + corners.col(1)) * 0.5;
//             // // evaluate a at midpoint
//             // auto a_eval_sqrt_n = a_coeff_(cell, midpoint)[0].sqrt() * normal;
//             // double A_F_sqrt = std::max(a_eval_sqrt_n[0], a_eval_sqrt_n[1]);
//             // double A_F = A_F_sqrt * A_F_sqrt;

//             //discontinuity penalization
//             auto disc_pen = disc_pen_(*edge, A_F);

            

//             //!!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!!
//             if (!boundary_edge_(*edge)){
                
//                 //pointers to adjacent polygons and the edge's local index in those
//                 auto polygon_plus = polygon_pair.first.first;
//                 auto polygon_minus = polygon_pair.second.first;
//                 auto edge_sub_idx_plus = polygon_pair.first.second;
//                 auto edge_sub_idx_minus = polygon_pair.second.second;

//                 //dof info
//                 nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
//                 nonstd::span<const Eigen::Index> dof_minus(dofhandler.GlobalDofIndices(*polygon_minus));

//                 //local - global mappings
//                 lf::dgfe::BoundingBox box_plus(*polygon_plus);
//                 lf::dgfe::BoundingBox box_minus(*polygon_minus);

//                 //mapped qr points to respective bounding box
//                 Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
//                 Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};


//                 //Assemble wi, vi, wj, vj
//                 Eigen::VectorXd wi = Eigen::VectorXd::Zero(n_basis);
//                 Eigen::VectorXd vi = Eigen::VectorXd::Zero(n_basis);
//                 Eigen::VectorXd wj = Eigen::VectorXd::Zero(n_basis);
//                 Eigen::VectorXd vj = Eigen::VectorXd::Zero(n_basis);
//                 //loop over basis functions
//                 for (int basis = 0; basis < n_basis; basis++){

//                         //loop over qr points
//                         for (int i = 0; i< gram_dets_s.size(); i++){
//                             //wi and vi
//                             double wi_p = legendre_basis(basis, max_legendre_degree_, zeta_box_plus.col(i)) * w_ref_s[i] * gram_dets_s[i];
//                             wi[basis] += wi_p;
//                             vi[basis] += wi_p;

//                             //wj and vj
//                             double wj_p = legendre_basis(basis, max_legendre_degree_, zeta_box_minus.col(i)) * w_ref_s[i] * gram_dets_s[i];
//                             wj[basis] += wj_p;
//                             vj[basis] += wj_p;
//                         }
    
//                 }

//                 //assembly double loop
//                 for (int l = 0; l < n_basis; l++) {
//                     for (int k = 0; k < n_basis; k++) {
//                         //part wivi
//                         matrix.AddToEntry(dof_plus[l], dof_plus[k], wi[l] * vi[k] * disc_pen);
//                         //part wivj
//                         matrix.AddToEntry(dof_plus[l], dof_minus[k], -wi[l] * vj[k] * disc_pen);
//                         //part wjvi
//                         matrix.AddToEntry(dof_minus[l], dof_plus[k], -wj[l] * vi[k] * disc_pen);
//                         //part wjvj
//                         matrix.AddToEntry(dof_minus[l], dof_minus[k], wj[l] * vj[k] * disc_pen);
//                     }
//                 }  // end assembly local double loop

//             }
//             //still second term but for boundary edges => simpler expression
//             if (boundary_d_edge_(*edge)){

//                 auto polygon_plus = polygon_pair.first.first;
//                 //dof info
//                 nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                
//                 //loop over basis functions in trial space
//                 for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
//                     //loop over bsis functions in test space
//                     for(int basis_test = 0; basis_test < n_basis; basis_test++){

//                         SCALAR sum = 0.0;
//                         //sum over qr points
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum * disc_pen);
//                     }
//                 }
//             }
//             //!!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!!


//             //!!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!!

//             //diffusion coefficient matrix evaluated at qr points
//             auto a = a_coeff_(cell, zeta_box_s);

//             if (!boundary_edge_(*edge)){
                
//                 //pointers to adjacent polygons and the edge's local index in those
//                 auto polygon_plus = polygon_pair.first.first;
//                 auto polygon_minus = polygon_pair.second.first;
//                 auto edge_sub_idx_plus = polygon_pair.first.second;
//                 auto edge_sub_idx_minus = polygon_pair.second.second;

//                 //dof info
//                 nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
//                 nonstd::span<const Eigen::Index> dof_minus(dofhandler.GlobalDofIndices(*polygon_minus));

//                 //local - global mappings
//                 lf::dgfe::BoundingBox box_plus(*polygon_plus);
//                 lf::dgfe::BoundingBox box_minus(*polygon_minus);

//                 //mapped qr points to respective bounding box
//                 Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
//                 Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};

//                 auto normal_plus = normal;
//                 auto normal_minus = -1.0 * normal_plus;


//                 //loop over basis functions in trial space
//                 for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
//                     //loop over bsis functions in test space
//                     for(int basis_test = 0; basis_test < n_basis; basis_test++){
                        
//                         //First, test functions on plus side are nonzero
//                         ///////////////////////////////////////////////
//                         SCALAR sum = 0.0;
    
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
//                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
//                                                               legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
//                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
//                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -0.5 * sum);

//                         //Secondly, test functions on minus side are nonzero
//                         ////////////////////////////////////////////////////
//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
//                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
//                                                               legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
//                                                              legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -0.5 * sum);

//                         sum = 0.0;
//                         for (int i = 0; i < gram_dets_s.size(); i++){
//                             Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
//                                                              legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

//                             sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_minus)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], -0.5 * sum);
//                     }
//                 }
//             }

//             //still third term but for boundary edges => simpler expression
//             if (boundary_d_edge_(*edge)){

//                 auto polygon_plus = polygon_pair.first.first;
//                 //dof info
//                 nonstd::span<const Eigen::Index> dof_plus(dofhandler.GlobalDofIndices(*polygon_plus));
                
//                 //loop over basis functions in trial space
//                 for (int basis_trial = 0; basis_trial < n_basis; basis_trial++){
//                     //loop over bsis functions in test space
//                     for(int basis_test = 0; basis_test < n_basis; basis_test++){

//                         SCALAR sum = 0.0;
//                         //sum over qr points
//                         for (int i = 0; i < gram_dets_s.size(); i++){
                            
//                             Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
//                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

//                             Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
//                                                             legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

//                             sum +=      (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * normal)
//                                     +   (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * normal)
//                                     * w_ref_s[i] * gram_dets_s[i];
//                         }
//                         matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -sum);
//                     }
//                 }
//             }
                


//     }


//     }

// private:

//     std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
//     unsigned integration_degree_;
//     unsigned max_legendre_degree_;
//     lf::quad::QuadRuleCache qr_cache_;
//     DIFFUSION_COEFF a_coeff_;
//     EDGESELECTOR boundary_edge_;
//     EDGESELECTOR boundary_d_edge_;
//     lf::dgfe::DiscontinuityPenalization disc_pen_;
//     //used to make sure the second and third term is only evaluated once per edge
//     //is true if edge has already been evaluated
//     EDGESELECTOR evaluated_edge_;
//     l2_proj_sqrt_a_nabla_basis &l2_projection_;

//     lf::mesh::utils::CodimMeshDataSet<bool> initialize_evaluated_edge(){
//         lf::mesh::utils::CodimMeshDataSet<bool> result(dgfe_space_ptr_->Mesh(), 1, false);
//         return result;
//     }
// };

// } //namespace lf::dgfe

// #endif // DIFFUSION_DGFE_H