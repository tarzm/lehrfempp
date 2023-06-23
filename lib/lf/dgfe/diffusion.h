
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

        //declare dof_plus and dof_minus for lambda capture
        nonstd::span<const Eigen::Index> dof_plus;
        nonstd::span<const Eigen::Index> dof_minus;
        //declare basis_test and basis_trial for lambda capture
        int basis_test;
        int basis_trial;

        //lambda for debugging
        auto galerkin_debug = [&dof_plus, &dof_minus, &basis_trial, &basis_test](int row, int col, double value, bool test_plus, bool trial_plus,
                                std::string additional = "", double additional_value = 0) -> void {
            
            int dof_row = test_plus ? dof_plus[basis_test] : dof_minus[basis_test];
            int dof_col = trial_plus ? dof_plus[basis_trial] : dof_minus[basis_trial];

            if (row == dof_row && col == dof_col && !(additional == "")){
                std::cout << "\tAdded " << value << " to G(" << row << " , " << col << ") with" << additional << " = \n\t\t\t" << additional_value << "\n";
            } else if(row == dof_row && col == dof_col){
                std::cout << "\tAdded " << value << " to G(" << row << " , " << col << ")\n";
            }
        };
        
        
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

            //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!

            auto cell_global_idx = dgfe_space_ptr_->Mesh()->Index(*cell);
            std::cout << "\nCell " << cell_global_idx << " in term 1\n";

            //dof info
            dof_plus = dofhandler.GlobalDofIndices(*cell);

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

                //diffusion coefficient matrix evaluated at qr points
                auto a = a_coeff_(*cell, zeta_box_t);

                //loop over basis functions in trial space
                for (basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over basis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){
                        //sum over qr points
                        SCALAR sum = 0;
                        for (int i = 0; i < gram_dets_t.size(); i++){

                            Eigen::Vector2d nabla_w{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                                    legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};
                            Eigen::Vector2d nabla_v{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(0),
                                                    legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * box.inverseJacobi(1)};

                            sum += (a[i] * nabla_w).dot(nabla_v) * w_ref_t[i] * gram_dets_t[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, sum, 1, 1);
                        //DEBUG
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum);
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
            //std::cout << "Zeta Global S:\n\t\t" << zeta_global_s << "\n";
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
            //gramian determinants
            Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};


            //CHECK ORIENTATION INTAGRATION
            //coordinates of nodes
            auto node_0 = edge->SubEntities(1)[0];
            auto node_1 = edge->SubEntities(1)[1];
            auto corners_0 = lf::geometry::Corners(*(node_0->Geometry()));
            auto corners_1 = lf::geometry::Corners(*(node_1->Geometry()));

            auto distance_0 = (zeta_global_s.col(0)-corners_0).norm();
            auto distance_1 = (zeta_global_s.col(0)-corners_1).norm();
            if (distance_0 < distance_1){
                std::cout << "RIGHT\n";
            } else {
                std::cout << "WRONG\n";
            }
            //CHECK ORIENTATION INTAGRATION

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

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                std::cout << "\nEdge " << edge_global_idx << " in term 2\n";
                
                //pointers to adjacent polygons and the edge's local index in those
                auto polygon_plus = polygon_pair.first.first;
                auto polygon_minus = polygon_pair.second.first;
                auto edge_sub_idx_plus = polygon_pair.first.second;
                auto edge_sub_idx_minus = polygon_pair.second.second;



                //dof info
                dof_plus = dofhandler.GlobalDofIndices(*polygon_plus);
                dof_minus = dofhandler.GlobalDofIndices(*polygon_minus);

                //local - global mappings
                lf::dgfe::BoundingBox box_plus(*polygon_plus);
                lf::dgfe::BoundingBox box_minus(*polygon_minus);

                //mapped qr points to respective bounding box
                Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
                Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};

                //loop over basis functions in trial space
                for (basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over basis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){
                        

                        SCALAR sum = 0.0;
                        //First part wi+ * vi+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, sum*disc_pen, 1, 1);
                        //DEBUG
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum * disc_pen);
                        
                        sum = 0.0;
                        //Second part - wi+ * vj+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, -sum*disc_pen, 0, 1);
                        //DEBUG
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], -sum * disc_pen);

                        sum = 0.0;
                        //Third part - wj+ * vi+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, -sum*disc_pen, 1, 0);
                        //DEBUG
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -sum * disc_pen);

                        sum = 0.0;
                        //Third part + wj+ * vj+
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, sum*disc_pen, 0, 0);
                        //DEBUG
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], sum * disc_pen);
                    }
                }
            }
            //still second term but for boundary edges => simpler expression
            if (boundary_d_edge_(*edge)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                std::cout << "\nEdge " << edge_global_idx << " in term 2\n";

                auto polygon_plus = polygon_pair.first.first;
                //dof info
                dof_plus = dofhandler.GlobalDofIndices(*polygon_plus);
                
                //loop over basis functions in trial space
                for (basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){

                        SCALAR sum = 0.0;
                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            sum += legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, sum*disc_pen, 1, 1);
                        //DEBUG
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], sum * disc_pen);
                    }
                }
            }
            //!!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!!
            

            //!!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!!

            //diffusion coefficient matrix evaluated at qr points

            auto a = a_coeff_(cell, zeta_box_s);


            if (!boundary_edge_(*edge)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                std::cout << "\nEdge " << edge_global_idx << " in term 3\n";


                
                //pointers to adjacent polygons and the edge's local index in those
                auto polygon_plus = polygon_pair.first.first;

                auto polygon_minus = polygon_pair.second.first;
                auto edge_sub_idx_plus = polygon_pair.first.second;
                auto edge_sub_idx_minus = polygon_pair.second.second;


                //dof info
                dof_plus = dofhandler.GlobalDofIndices(*polygon_plus);
                dof_minus = dofhandler.GlobalDofIndices(*polygon_minus);


                //local - global mappings
                lf::dgfe::BoundingBox box_plus(*polygon_plus);
                lf::dgfe::BoundingBox box_minus(*polygon_minus);

                //mapped qr points to respective bounding box
                Eigen::MatrixXd zeta_box_plus{box_plus.inverseMap(zeta_global_s)};
                Eigen::MatrixXd zeta_box_minus{box_minus.inverseMap(zeta_global_s)};

                auto normal_plus = normal;
                auto normal_minus = -1.0 * normal_plus;

                //loop over basis functions in trial space
                for (basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){
                        
                        //First part {{a * nabla_w}} * [[v]]
                        ///////////////////////////////////////////////
                        SCALAR sum = 0.0;



                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, -0.5*sum, 1, 1);
                        //DEBUG
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -0.5 * sum);
 
                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, 0.5*sum, 0, 1);
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial], 0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, -0.5*sum, 1, 0, " Here 1", 1);
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_trial_minus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_trial_minus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        //DEBUG
                        galerkin_debug(18, 30, 0.5*sum, 0, 0);
                        //DEBUG
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], 0.5 * sum);

                        //Second part {{a * v}} * [[w]]
                        ////////////////////////////////////////////////////
                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, -0.5*sum, 1, 1);
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], - 0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_plus.col(i)) * box_plus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, 0.5*sum, 0, 1, "Here 2", 2);
                        matrix.AddToEntry(dof_plus[basis_test], dof_minus[basis_trial], 0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_plus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, -0.5*sum, 1, 0);
                        matrix.AddToEntry(dof_minus[basis_test], dof_plus[basis_trial],  -0.5 * sum);

                        sum = 0.0;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            Eigen::Vector2d nabla_test_minus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(0),
                                                              legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_minus.col(i)) * box_minus.inverseJacobi(1)};

                            sum += (a[i] * nabla_test_minus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_minus.col(i)) * normal_plus)
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, 0.5*sum, 0, 0);
                        matrix.AddToEntry(dof_minus[basis_test], dof_minus[basis_trial], 0.5 * sum);
                    }
                }
            }

            //still third term but for boundary edges => simpler expression
            if (boundary_d_edge_(*edge)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                std::cout << "\nEdge " << edge_global_idx << " in term 3\n";

                auto polygon_plus = polygon_pair.first.first;
                //dof info
                dof_plus = dofhandler.GlobalDofIndices(*polygon_plus);
                
                //loop over basis functions in trial space
                for (basis_trial = 0; basis_trial < n_basis; basis_trial++){
                    //loop over bsis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){

                        SCALAR sum = 0.0;
                        //sum over qr points
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            
                            Eigen::Vector2d nabla_trial_plus{legendre_basis_dx(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                             legendre_basis_dy(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            // if (basis_trial == 2 && basis_test == 0){
                            //     std::cout << "Nabla_basis_trial = \n" << nabla_trial_plus << "\n";
                            // }

                            Eigen::Vector2d nabla_test_plus{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                            legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                            // if (basis_trial == 2 && basis_test == 0){
                            //     std::cout << "Nabla_basis_test = \n" << nabla_test_plus << "\n";

                            //     std::cout << "Normal = \n" << normal << "\n";

                            //     std::cout << "basis_test = \n" << legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) << "\n";

                            //     std::cout << "A = \n" << a[i] << "\n";

                            //     std::cout << "Erster Teil = \n" << (a[i] * nabla_trial_plus).dot(
                            //                         legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * normal)* w_ref_s[i] * gram_dets_s[i] << "\n";

                            //     std::cout << "Zweiter Teil = \n" << (a[i] * nabla_test_plus).dot(
                            //                         legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * normal)* w_ref_s[i] * gram_dets_s[i] << "\n";
                            // }   

                            sum +=      ((a[i] * nabla_trial_plus).dot(legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * normal)
                                    +   (a[i] * nabla_test_plus).dot(legendre_basis(basis_trial, max_legendre_degree_, zeta_box_s.col(i)) * normal))
                                    * w_ref_s[i] * gram_dets_s[i];
                        }
                        galerkin_debug(18, 30, -sum, 1, 1);
                        matrix.AddToEntry(dof_plus[basis_test], dof_plus[basis_trial], -sum);
                    }
                }
            }
            
            //!!!!!!!!!!!!! END THIRD TERM !!!!!!!!!!!!!!
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