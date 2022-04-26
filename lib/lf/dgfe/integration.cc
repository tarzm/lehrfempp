/**
 * @file
 * @brief Implementation of functions for integrating polynomials on polygons in 2D
 * @author Tarzis Maurer
 * @date 2022-19-04
 * @copyright ETH Zurich
 */

#include "integration.h"


namespace lf::dgfe {


//it is assumed that all edges are oriented counterclockwise
//such that edge.col(1) is the coordinate of the second point
Eigen::MatrixXd outwardNormal(const Eigen::MatrixXd edge){
    LF_VERIFY_MSG( edge.cols() == 2 && edge.rows() == 2, "Edge whose outward normal should be found must be of shape 2x2");
    Eigen::Vector3d tangent(edge(0,1) - edge(0,0), edge(1,1) - edge(1,0), 0.0);
    Eigen::Vector3d z(0.0, 0.0, 1);
    return (tangent.cross(z)).head<2>().normalized();
}

scalar_t euclideanDist(const Eigen::MatrixXd a, const Eigen::MatrixXd b){
    LF_VERIFY_MSG( a.cols() == 1 && a.rows() == 2, "Matrices (vectors) must be of shape 2x1");
    LF_VERIFY_MSG( b.cols() == 1 && b.rows() == 2, "Matrices (vectors) must be of shape 2x1");
    Eigen::MatrixXd a_b(2,1);
    a_b <<  b(0,0) - a(0,0), 
            b(1,0) - a(1,0);
    return a_b.norm();
}

scalar_t integrate(Eigen::MatrixXd corners, int degree_x, int degree_y){

    auto n_cols = corners.cols();
    switch(n_cols){
        case 1: {       // entity is a point
            //std::cout << "POINT: " << std::pow(corners(0,0), degree_x) * std::pow(corners(1,0), degree_y) << "\n";
            return std::pow(corners(0,0), degree_x) * std::pow(corners(1,0), degree_y);
        }

        case 2: {       // entity is an edge
            scalar_t sum = 0.0;
            scalar_t pre_factor = 1.0 / (1.0 + degree_x + degree_y);

            //setup of x_0
            scalar_t x0_r;
            size_type r = (degree_x < degree_y) ? 0 : 1;
            auto n = outwardNormal(corners);
            scalar_t b = (n.col(0)).dot(corners.col(0));
            //r cannot be the axis where n_r is 0 => if so, flip it
            if (n(r,0) == 0.0){
                r = 1 - r;
            }
            x0_r = b / n(r,0);
            Eigen::MatrixXd x0 = Eigen::MatrixXd::Zero(2,1);
            x0(r,0) = x0_r;

            //run over edge's end points
            for (int i = 0; i < 2; i++){
                //term in sum
                sum += euclideanDist(corners.col(i), x0) * integrate(corners.col(i), degree_x, degree_y);
            }

            //term x0_i * k_i * I(N, ...)
            if (r == 0 && degree_x > 0){
                sum += x0(r,0) * degree_x * integrate(corners, degree_x - 1, degree_y);
            } else if (r == 1 && degree_y > 0) {
                sum += x0(r,0) * degree_y * integrate(corners, degree_x, degree_y - 1);
            }

            sum *= pre_factor;
            //std::cout << "EDGE: " << sum << "\n";
            return sum;
        }

        default: {              //entitiy is a polygon
            scalar_t sum = 0.0;
            scalar_t pre_factor = 1.0 / (2.0 + degree_x + degree_y);

            for (int i = 0; i < corners.cols(); i++){
                Eigen::MatrixXd edge(2,2);
                edge.col(0) = corners.col(i);
                edge.col(1) = corners.col((i + 1) % corners.cols());
                scalar_t b_i =  (edge.col(0))(0,0) * (outwardNormal(edge))(0,0) +
                                (edge.col(0))(1,0) * (outwardNormal(edge))(1,0); //basically dot product
                sum += b_i * integrate(edge, degree_x, degree_y);
            }
            sum *= pre_factor;
            //std::cout << "POLYGON: " << sum << "\n";
            return sum;
        }
    }
}


} //namespace lf::dgfe