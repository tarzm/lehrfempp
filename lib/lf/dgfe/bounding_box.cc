/**
 * @file
 * @brief Functionalities of the axis aligned bounding box
 * @author Tarzis Maurer
 * @date 2022-03-05
 * @copyright ETH Zurich
 */

#include "bounding_box.h"


namespace lf::dgfe {

BoundingBox::BoundingBox(const lf::mesh::Entity &entity){
    Eigen::MatrixXd corners;
    //retrieve coordniates of vertices of entity
    switch(entity.RefEl()){
        case lf::base::RefEl::kPoint():
        case lf::base::RefEl::kSegment():
            LF_VERIFY_MSG(true, "Entity must be of codim 0");
            break;
        
        case lf::base::RefEl::kTria():
        case lf::base::RefEl::kQuad():
            corners = lf::geometry::Corners(*entity.Geometry());
            break;

        case lf::base::RefEl::kPolygon():
            corners = lf::mesh::polytopic2d::Corners(&entity);
            break;
    }
    LF_VERIFY_MSG(corners.rows() == 2, "Coordinates must be of 2 dimensions");

    Eigen::MatrixXd bounding(2,2);
    bounding.col(0) = corners.rowwise().minCoeff();
    bounding.col(1) = corners.rowwise().maxCoeff();

    translation_ = bounding.rowwise().mean();
    jacobi_ = Eigen::Matrix2d::Zero();
    jacobi_(0,0) = (bounding(0,1) - bounding(0,0)) / 2.0;
    jacobi_(1,1) = (bounding(1,1) - bounding(1,0)) / 2.0;

    inverse_jacobi_ = jacobi_.inverse();
}

Eigen::MatrixXd BoundingBox::map(const Eigen::MatrixXd corners){
    LF_VERIFY_MSG(corners.rows() == 2, "Coordinates must be of 2 dimensions");
    Eigen::MatrixXd translation_matrix(corners.rows(), corners.cols());
    for (int i = 0; i < corners.cols(); i++){
        translation_matrix.col(i) = translation_;
    }
    return jacobi_ * corners + translation_matrix;
}

Eigen::MatrixXd BoundingBox::inverseMap(const Eigen::MatrixXd corners){
    LF_VERIFY_MSG(corners.rows() == 2, "Coordinates must be of 2 dimensions");
    Eigen::MatrixXd translation_matrix(corners.rows(), corners.cols());
    for (int i = 0; i < corners.cols(); i++){
        translation_matrix.col(i) = translation_;
    }
    return inverse_jacobi_ * (corners - translation_matrix);
}

scalar_t BoundingBox::det(){
    return jacobi_(0,0) * jacobi_(1,1);
}

scalar_t BoundingBox::inverseJacobi(unsigned i){
    LF_VERIFY_MSG(i == 0 || i == 1, "Requested entry must be (0, 0) or (1, 1)");
    return inverse_jacobi_(i,i);
}

} //namespace lf::dgfe