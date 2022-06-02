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



namespace lf::dgfe {

/**Type of scalar used*/
using scalar_t = double;

using size_type = lf::mesh::polytopic2d::Mesh::size_type;

//it is assumed that all edges are oriented counterclockwise
//such that edge.col(1) is the coordinate of the second point
Eigen::MatrixXd outwardNormal(const Eigen::MatrixXd edge);

scalar_t euclideanDist(const Eigen::MatrixXd a, const Eigen::MatrixXd b);

scalar_t integrate(Eigen::MatrixXd corners, int degree_x, int degree_y);

} //namespace dgfe

#endif // LF_DGFEINTEGRATION_H