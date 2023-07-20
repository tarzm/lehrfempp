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

#include <unsupported/Eigen/MatrixFunctions>

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"
#include "integration.h"
#include "mesh_function_dgfe.h"
#include "mesh_function_global.h"
#include "discontinuity_penalization.h"

namespace lf::dgfe {

template <typename Derived>
std::string get_shape(const Eigen::EigenBase<Derived>& x)
{
    std::ostringstream oss;
    oss  << "(" << x.rows() << ", " << x.cols() << ")";
    return oss.str();
}


template<typename SCALAR, typename DIFFUSION_COEFF, typename ADVECTION_COEFF, typename EDGESELECTOR, typename MESHFUNC_F, typename MESHFUNC_gD, typename MESHFUNC_gN, typename TMPMATRIX>
class AdvectionReactionDiffusionRHSAssembler {

public:

    using l2_proj_sqrt_a_nabla_basis = std::pair<std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>, std::vector<lf::dgfe::MeshFunctionDGFE<SCALAR>>>;

    AdvectionReactionDiffusionRHSAssembler(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, MESHFUNC_F f, MESHFUNC_gD gD, MESHFUNC_gN gN,
                                    DIFFUSION_COEFF a_coeff, ADVECTION_COEFF b_coeff, 
                                    EDGESELECTOR boundary_minus_edge, EDGESELECTOR boundary_d_edge,
                                    EDGESELECTOR boundary_n_edge, unsigned integration_degree, lf::dgfe::DiscontinuityPenalization disc_pen,
                                    l2_proj_sqrt_a_nabla_basis &l2_proj)
        : dgfe_space_ptr_(std::move(dgfe_space_ptr)), integration_degree_(integration_degree),
         max_legendre_degree_(dgfe_space_ptr_->MaxLegendreDegree()), b_coeff_(b_coeff), a_coeff_(a_coeff),
         boundary_minus_edge_(std::move(boundary_minus_edge)), boundary_d_edge_(std::move(boundary_d_edge)),
         boundary_n_edge_(std::move(boundary_n_edge)), f_(f), gD_(gD), gN_(gN), disc_pen_(disc_pen), l2_projection_(l2_proj) {
            LF_VERIFY_MSG(dgfe_space_ptr_ != nullptr, "No DGFE space defined");
            LF_VERIFY_MSG(dgfe_space_ptr_ == disc_pen_.dgfeSpace(), "Space in constructor and space of discontinuity penalization do not match");
    }


    void assemble(TMPMATRIX &rhs_vec){
        
        nonstd::span<const Eigen::Index> dofs_cell;
        int basis_test;

        //lambda for debugging
        auto rhs_debug = [&dofs_cell, &basis_test](int row, double value, std::string additional = "", double additional_value = 0.0) -> void {
            
            int dof_row = dofs_cell[basis_test];

            if (row == dof_row && !(additional == "")){
                // std::cout << "\tAdded " << value << " to RHS " << row << " with " << additional << " = \n\t\t\t" << additional_value << "\n";
            } else if(row == dof_row){
                // std::cout << "\tAdded " << value << " to RHS " << row << "\n";
            }
        };
        

        //number of basis functions per cell
        const unsigned n_basis = (max_legendre_degree_ == 1) ? 4 : 9;
        //alias for dofhandler
        auto dofhandler = dgfe_space_ptr_->LocGlobMap();

        //////////quadrule setup triangle
        const lf::quad::QuadRule qr_t = qr_cache_.Get(lf::base::RefEl::kTria(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};

        ////////////quadrule setup segment
        const lf::quad::QuadRule qr_s = qr_cache_.Get(lf::base::RefEl::kSegment(), integration_degree_);
        // qr points
        const Eigen::MatrixXd zeta_ref_s{qr_s.Points()};
        //weights
        Eigen::VectorXd w_ref_s{qr_s.Weights()};



        //!!!!!!!!!!!!! FIRST TERM !!!!!!!!!!!!!!
        //     f * v   over all cells

        //loop over cells
        for (const lf::mesh::Entity *cell : dgfe_space_ptr_->Mesh()->Entities(0)){

            auto cell_global_idx = dgfe_space_ptr_->Mesh()->Index(*cell);
                //std::cout << "\nCell " << cell_global_idx << " in term 1\n";
        

            //local - global mapping
            lf::dgfe::BoundingBox box(*cell);
            //get sub-tessellation
            auto sub_tessellation = subTessellation(cell);
            //dof info
            dofs_cell = dofhandler.GlobalDofIndices(*cell);

            
            //loop over triangles in the sub-tessellation
            for(auto& tria_geo_ptr : sub_tessellation){
                // qr points mapped to triangle
                Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
                //gramian determinants
                Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

                auto f_evaluated = f_(*cell, zeta_box_t);

                scalar_t sum;
                //loop over basis functions in test space
                for (basis_test = 0; basis_test < n_basis; basis_test++){
                    sum = 0;
                    //sum over qr points
                    for (int i = 0; i < gram_dets_t.size(); i++){
                        
                        sum += f_evaluated[i] * legendre_basis(basis_test, max_legendre_degree_, zeta_box_t.col(i)) * w_ref_t[i] * gram_dets_t[i];

                    }
                    rhs_debug(1,sum);
                    rhs_vec[dofs_cell[basis_test]] += sum;
                }
            }
            //!!!!!!!!!!!!! END FIRST TERM !!!!!!!!!!!!!

            //!!!!!!!!!!!!! SECOND TERM !!!!!!!!!!!!!

            size_type edge_sub_idx = 0;

            //loop over cell's edges
            for (auto edge : cell->SubEntities(1)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                //std::cout << "\nEdge " << edge_global_idx << " in term 2\n";

                //normal n
                auto normal = lf::dgfe::outwardNormal(lf::geometry::Corners(*(edge->Geometry())));
                //if orientation of edge in polygon is negative, normal has to be multiplied by -1;
                normal *= (int) (cell->RelativeOrientations()[edge_sub_idx]);

                // qr points mapped to segment
                Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
                //gramian determinants
                Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};

                auto b = b_coeff_(*cell, zeta_box_s);

                //does edge belong do delta_minus_kappa
                SCALAR result_delta_minus_kappa = 0.0;
                for (int i = 0; i < gram_dets_s.size(); i++){
                    result_delta_minus_kappa += b[i].dot(normal) * w_ref_s[i] * gram_dets_s[i];
                }
                bool delta_minus_kappa = result_delta_minus_kappa < 0 ? true : false; 

                //    - ( b * n ) * gD * v^+   over all edges which are either on boundary_D or boundary_minus and belong to delta_minus_k of the cell
                if (delta_minus_kappa && (boundary_d_edge_(*edge) || boundary_minus_edge_(*edge))){

                    auto gD_evaluated = gD_(*cell, zeta_box_s);
                    
                    //loop over bsis functions in test space
                    for(basis_test = 0; basis_test < n_basis; basis_test++){
                        //sum over qr points and add if delta_minus_k
                        scalar_t sum;
                        for (int i = 0; i < gram_dets_s.size(); i++){
                            if (b[i].dot(normal) < 0){
                                sum -= (b[i].dot(normal)) * gD_evaluated[i]
                                                                    * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                                                    * w_ref_s[i] * gram_dets_s[i];
                            }
                        }
                        rhs_debug(1,sum);
                        rhs_vec[dofs_cell[basis_test]] += sum;
                    }
                }
                edge_sub_idx++;
            }
        }
        //!!!!!!!!!!!!! END SECOND TERM !!!!!!!!!!!!!!              
        
        //loop over edges
        for (const lf::mesh::Entity *edge : dgfe_space_ptr_->Mesh()->Entities(1)){

            //!!!!!!!!!!!!! THIRD TERM !!!!!!!!!!!!!!
            //    -(a * nabla(v) * n - disc_pen * v) dS over boundary_d
            if (boundary_d_edge_(*edge)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                //std::cout << "\nEdge " << edge_global_idx << " in term 3\n";

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

                //dof info
                dofs_cell = dofhandler.GlobalDofIndices(cell);
                
                //calculate A_F
                Eigen::MatrixXd A_F_mat = Eigen::MatrixXd::Zero(2, gram_dets_s.size());
                for (int i = 0; i < gram_dets_s.size(); i++){
                    A_F_mat.col(i) = a_coeff_(cell, zeta_box_s.col(i))[0].sqrt() * normal;
                }
                SCALAR A_F = A_F_mat.lpNorm<Eigen::Infinity>();
                A_F = A_F * A_F;

                //discontinuity penalty
                scalar_t disc_penalty = disc_pen_(*edge, A_F);

                //gD
                auto gD_evaluated = gD_(cell, zeta_box_s);

                //loop over basis functions in test space
                for(basis_test = 0; basis_test < n_basis; basis_test++){
                    //sum over qr points
                    scalar_t sum = 0;
                    for (int i = 0; i < gram_dets_s.size(); i++){
                        Eigen::Vector2d nabla_v{legendre_basis_dx(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(0),
                                                legendre_basis_dy(basis_test, max_legendre_degree_, zeta_box_s.col(i)) * box.inverseJacobi(1)};

                        sum -= gD_evaluated[i] * ( (a_coeff_(cell, zeta_box_s.col(i))[0] *  nabla_v).dot(normal)
                                                            - disc_penalty * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i)) )
                                                * w_ref_s[i] * gram_dets_s[i];
                    }
                    rhs_debug(1,sum);
                    rhs_vec[dofs_cell[basis_test]] += sum;
                }          
            }
            //!!!!!!!!!!!!! END THIRD TERM !!!!!!!!!!!!!!


            //!!!!!!!!!!!!! FOURTH TERM !!!!!!!!!!!!!!
            //    + gN * v dS over boundary_N
            if (boundary_n_edge_(*edge)){

                auto edge_global_idx = dgfe_space_ptr_->Mesh()->Index(*edge);
                //std::cout << "\nEdge " << edge_global_idx << " in term 4\n";

                //adjacent cell
                auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(edge);
                const lf::mesh::Entity &cell = *polygon_pair.first.first;
                
                //bounding box
                lf::dgfe::BoundingBox box(cell);
                // qr points mapped to segment
                Eigen::MatrixXd zeta_global_s{edge->Geometry()->Global(zeta_ref_s)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_s{box.inverseMap(zeta_global_s)};
                //gramian determinants
                Eigen::VectorXd gram_dets_s{edge->Geometry()->IntegrationElement(zeta_ref_s)};

                //dof info
                nonstd::span<const Eigen::Index> dofs_cell(dofhandler.GlobalDofIndices(cell));
                
                //gD
                auto gN_evaluated = gN_(cell, zeta_box_s);

                //loop over basis functions in test space
                for(basis_test = 0; basis_test < n_basis; basis_test++){
                    //sum over qr points
                    scalar_t sum = 0;
                    for (int i = 0; i < gram_dets_s.size(); i++){
                        sum += gN_evaluated[i] * legendre_basis(basis_test, max_legendre_degree_, zeta_box_s.col(i))
                                                            * w_ref_s[i] * gram_dets_s[i];
                    }
                    rhs_debug(1,sum);
                    rhs_vec[dofs_cell[basis_test]] += sum;
                }
            }
            //!!!!!!!!!!!!! END FOURTH TERM !!!!!!!!!!!!!!
        }
    }

private:
    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
    unsigned integration_degree_;
    unsigned max_legendre_degree_;
    lf::quad::QuadRuleCache qr_cache_;
    MESHFUNC_F f_;
    MESHFUNC_gD gD_;
    MESHFUNC_gN gN_;
    DIFFUSION_COEFF a_coeff_;
    ADVECTION_COEFF b_coeff_;
    EDGESELECTOR boundary_minus_edge_;
    EDGESELECTOR boundary_d_edge_;
    EDGESELECTOR boundary_n_edge_;
    lf::dgfe::DiscontinuityPenalization disc_pen_;
    l2_proj_sqrt_a_nabla_basis &l2_projection_;
};



} //namespace lf::dgfe

#endif //ADVECTION_REACTION_DIFFUSION_RHS_H