#ifndef LF_DGFEINTEGRATION_H
#define LF_DGFEINTEGRATION_H

/**
 * @file
 * @brief Functions for integrating polynomials on polygons in 2D
 * @author Tarzis Maurer
 * @date 2022-19-04
 * @copyright ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/quad/quad.h>

#include "sub_tessellation.h"

namespace lf::dgfe {

/**Type of scalar used*/
using scalar_t = double;

using size_type = lf::mesh::polytopic2d::Mesh::size_type;

//it is assumed that all edges are oriented counterclockwise
//such that edge.col(1) is the coordinate of the second point
Eigen::MatrixXd outwardNormal(const Eigen::MatrixXd edge);

scalar_t euclideanDist(const Eigen::MatrixXd a, const Eigen::MatrixXd b);

scalar_t integrate(Eigen::MatrixXd corners, int degree_x, int degree_y);



template<typename SCALAR, typename FUNCTOR>
class SubTessellationIntegrator{
public:

    SubTessellationIntegrator(){
        for (int i = 0; i < 6; i++){
            segment_quadrules_[i] = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), i);
            tria_quadrules_[i] = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), i);
        }
    }


    SCALAR integrate(const lf::mesh::Entity *entity, FUNCTOR f, unsigned max_degree){
        auto ref_el = entity->RefEl();
        SCALAR sum = 0.0;
        lf::quad::QuadRule qr;

        switch(ref_el){
            case lf::base::RefEl::kPoint():{
                return f(lf::geometry::Corners(*(entity->Geometry())));
            }

            case lf::base::RefEl::kSegment():{
                //setup quadrature rule
                if (max_degree < 6){
                    qr = segment_quadrules_[max_degree];
                } else {
                    qr = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), max_degree);
                }
                //get geometry pointer
                auto edge_geo_ptr = entity->Geometry();
                // qr points
                const Eigen::MatrixXd zeta_ref{qr.Points()};
                // qr points mapped to segment
                Eigen::MatrixXd zeta{edge_geo_ptr->Global(zeta_ref)};
                //weights
                Eigen::VectorXd w_ref{qr.Weights()};
                //gramian determinants
                Eigen::VectorXd gram_dets{edge_geo_ptr->IntegrationElement(zeta_ref)};

                //sum over qr points
                for (int i = 0; i < qr.Points().cols(); i++){
                    sum += w_ref[i] * f(zeta.col(i)) * gram_dets[i];
                }
                return sum;
            }

            case lf::base::RefEl::kPolygon():{
                //setup quadrature rule
                if (max_degree < 6){
                    qr = tria_quadrules_[max_degree];
                } else {
                    qr = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), max_degree);
                }
                //get sub-tessellation
                auto sub_tessellation = subTessellation(entity);
                // qr points
                const Eigen::MatrixXd zeta_ref{qr.Points()};
                //weights
                Eigen::VectorXd w_ref{qr.Weights()};

                //loop over triangles in the sub-tessellation
                for (auto& tria_geo_ptr : sub_tessellation){
                    // qr points mapped to triangle
                    Eigen::MatrixXd zeta{tria_geo_ptr->Global(zeta_ref)};
                    //gramian determinants
                    Eigen::VectorXd gram_dets{tria_geo_ptr->IntegrationElement(zeta_ref)};
                    //sum over qr points
                    for (int i = 0; i < qr.Points().cols(); i++){
                        sum += w_ref[i] * f(zeta.col(i)) * gram_dets[i];
                    }
                }
                return sum;
            }

            default:{
                LF_ASSERT_MSG(false, "SubTessellation integration only implemented in the polytopic setting");
                return 0;
            }
        }
    }
private:
    std::array<lf::quad::QuadRule, 6> segment_quadrules_{};
    std::array<lf::quad::QuadRule, 6> tria_quadrules_{};

};








} //namespace dgfe

#endif // LF_DGFEINTEGRATION_H