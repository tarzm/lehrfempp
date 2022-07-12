/**
 * @file
 * @brief Implementation of functions for integrating polynomials on polygons in 2D
 * @author Tarzis Maurer
 * @date 2022-19-04
 * @copyright ETH Zurich
 */

#include "integration.h"

#include <lf/geometry/geometry.h>


namespace lf::dgfe {


//it is assumed that all edges are oriented counterclockwise
//such that edge.col(1) is the coordinate of the second point
Eigen::VectorXd outwardNormal(const Eigen::MatrixXd edge){
    LF_VERIFY_MSG( edge.cols() == 2 && edge.rows() == 2, "Edge whose outward normal should be found must be of shape 2x2");
    Eigen::Vector3d tangent(edge(0,1) - edge(0,0), edge(1,1) - edge(1,0), 0.0);
    Eigen::Vector3d z(0.0, 0.0, 1);
    auto vec2 = (tangent.cross(z)).head<2>();
    return vec2.normalized();
}

//calculates euclidean distance between two coordinates
scalar_t euclideanDist(const Eigen::MatrixXd a, const Eigen::MatrixXd b){
    LF_VERIFY_MSG( a.cols() == 1 && a.rows() == 2, "Matrices (vectors) must be of shape 2x1");
    LF_VERIFY_MSG( b.cols() == 1 && b.rows() == 2, "Matrices (vectors) must be of shape 2x1");
    Eigen::MatrixXd a_b(2,1);
    a_b <<  b(0,0) - a(0,0), 
            b(1,0) - a(1,0);
    return a_b.norm();
}

scalar_t integrate(Eigen::MatrixXd &corners, int degree_x, int degree_y){

    auto n_cols = corners.cols();
    switch(n_cols){
        case 1: {       // entity is a point
            // std::cout << "POINT: " << std::pow(corners(0,0), degree_x) * std::pow(corners(1,0), degree_y) << "\n";
            return std::pow(corners(0,0), degree_x) * std::pow(corners(1,0), degree_y);
        }

        case 2: {       // entity is an edge
            scalar_t sum = 0.0;
            
            //############ QUADRATURE-BASED APPROACH
            //get pointer to segments geometry
            auto edge_geo_ptr = std::make_unique<lf::geometry::SegmentO1>(corners);

            auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), degree_x + degree_y + 1);

            // qr points
            const Eigen::MatrixXd zeta_ref{qr.Points()};
            // qr points mapped to segment
            Eigen::MatrixXd zeta{edge_geo_ptr->Global(zeta_ref)};
            //weights
            Eigen::VectorXd w_ref{qr.Weights()};
            //gramian determinants
            Eigen::VectorXd gram_dets{edge_geo_ptr->IntegrationElement(zeta_ref)};

            for(int i = 0; i < qr.Points().cols(); i++){
                sum += w_ref[i] * (std::pow(zeta(0,i), degree_x) * std::pow(zeta(1,i), degree_y)) * gram_dets[i];
            }
            //############ QUADRATURE-BASED APPROACH END


            //############ ORIGINAL APPROACH ALGORITHM 1 FROM PAPER
            // scalar_t pre_factor = 1.0 / (1.0 + degree_x + degree_y);
            // Eigen::MatrixXd n = outwardNormal(corners);
            // scalar_t b = (n.col(0)).dot(corners.col(0));

            // //setup of x_0
            // scalar_t x0_r;
            // //r is where the exponent of the polynomial is smaller
            // size_type r = (degree_x <= degree_y) ? 0 : 1;
            // //r cannot be the axis where n_r is 0 => if so, flip it
            // if (n(r,0) == 0){
            //     r = 1 - r;
            // }
            // x0_r = b / n(r,0);
            // Eigen::MatrixXd x0 = Eigen::MatrixXd::Zero(2,1);
            // x0(r,0) = x0_r;

            // //run over edge's end points
            // for (int i = 0; i < 2; i++){
            //     //term in sum
            //     sum += euclideanDist(corners.col(i), x0) * integrate(corners.col(i), degree_x, degree_y);
            // }

            // //term x0_i * k_i * I(N, ...)
            // if (r == 0 && degree_x > 0){
            //     sum += x0(r,0) * degree_x * integrate(corners, degree_x - 1, degree_y);
            // } else if (r == 1 && degree_y > 0) {
            //     sum += x0(r,0) * degree_y * integrate(corners, degree_x, degree_y - 1);
            // }

            // sum *= pre_factor;
            //############ ORIGINAL APPROACH END

            
            // std::cout << "EDGE: " << sum << "\n";
            return sum;
        }

        default: {              //entitiy is a polygon
            scalar_t sum = 0.0;
            scalar_t pre_factor = 1.0 / (2.0 + degree_x + degree_y);

            for (int i = 0; i < corners.cols(); i++){
                Eigen::MatrixXd edge(2,2);
                edge.col(0) = corners.col(i);
                edge.col(1) = corners.col((i + 1) % corners.cols());
                scalar_t b_i = outwardNormal(edge).col(0).dot(edge.col(0));
                sum += b_i * integrate(edge, degree_x, degree_y);
            }
            sum *= pre_factor;
            // std::cout << "POLYGON: " << sum << "\n";
            return sum;
        }
    }
}


} //namespace lf::dgfe