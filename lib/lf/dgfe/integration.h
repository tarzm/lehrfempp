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
#include "bounding_box.h"

namespace lf::dgfe {

/**Type of scalar used*/
using scalar_t = double;

using size_type = lf::mesh::polytopic2d::Mesh::size_type;

//it is assumed that all edges are oriented counterclockwise
//such that edge.col(1) is the coordinate of the second point
Eigen::VectorXd outwardNormal(const Eigen::MatrixXd edge);

scalar_t euclideanDist(const Eigen::MatrixXd a, const Eigen::MatrixXd b);

scalar_t integrate(Eigen::MatrixXd corners, int degree_x, int degree_y);


template<typename SCALAR, typename MESHFUNC>
class SubTessellationIntegrator{
public:

    SubTessellationIntegrator(){};

    SCALAR integrate(const lf::mesh::Entity& entity, MESHFUNC f, unsigned max_degree){
        auto ref_el = entity.RefEl();
        SCALAR sum = 0.0;
        lf::quad::QuadRule qr;

        switch(ref_el){
            // case lf::base::RefEl::kPoint():{
            //     return f(*entity, lf::geometry::Corners(*(entity->Geometry())))[0];
            // }

            // case lf::base::RefEl::kSegment():{
            //     //setup quadrature rule
            //     qr = qr_cache_.Get(lf::base::RefEl::kSegment(), max_degree);
            //     //get geometry pointer
            //     auto edge_geo_ptr = entity->Geometry();
            //     // qr points
            //     const Eigen::MatrixXd zeta_ref{qr.Points()};
            //     // qr points mapped to segment
            //     Eigen::MatrixXd zeta{edge_geo_ptr->Global(zeta_ref)};
            //     //weights
            //     Eigen::VectorXd w_ref{qr.Weights()};
            //     //gramian determinants
            //     Eigen::VectorXd gram_dets{edge_geo_ptr->IntegrationElement(zeta_ref)};

            //     //sum over qr points
            //     for (int i = 0; i < qr.Points().cols(); i++){
            //         sum += w_ref[i] * f(entity, zeta.col(i)) * gram_dets[i];
            //     }
            //     return sum;
            // }

            case lf::base::RefEl::kPolygon():{
                //setup quadrature rule
                qr = qr_cache_.Get(lf::base::RefEl::kTria(), max_degree);
                //get sub-tessellation
                auto sub_tessellation = subTessellation(&entity);
                // qr points
                const Eigen::MatrixXd zeta_ref{qr.Points()};
                //weights
                Eigen::VectorXd w_ref{qr.Weights()};
                //bounding box
                lf::dgfe::BoundingBox box(entity);

                //loop over triangles in the sub-tessellation
                for (auto& tria_geo_ptr : sub_tessellation){
                    // qr points mapped to triangle
                    Eigen::MatrixXd zeta_global{tria_geo_ptr->Global(zeta_ref)};
                    // qr points mapped back to reference bounding box
                    Eigen::MatrixXd zeta_box{box.inverseMap(zeta_global)};
                    //gramian determinants
                    Eigen::VectorXd gram_dets{tria_geo_ptr->IntegrationElement(zeta_ref)};
                    
                    //function evaluated at the local points (points within the reference bounding box)
                    auto f_res_vec = f(entity, zeta_box);

                    Eigen::VectorXd f_result = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(f_res_vec.data(), f_res_vec.size());

                    sum += w_ref.cwiseProduct(f_result).dot(gram_dets);
                }
                return sum;
            }

            default:{
                LF_ASSERT_MSG(false, "SubTessellation integration only implemented for polygons");
                return 0;
            }
        }
    }
private:
    lf::quad::QuadRuleCache qr_cache_;

};

} //namespace dgfe

#endif // LF_DGFEINTEGRATION_H