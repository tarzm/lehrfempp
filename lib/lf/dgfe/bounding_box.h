#ifndef LF_BOUNDING_BOX_H
#define LF_BOUNDING_BOX_H

/**
 * @file
 * @brief Functionalities of the axis aligned bounding box
 * @author Tarzis Maurer
 * @date 2022-03-05
 * @copyright ETH Zurich
 */


#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <Eigen/Eigen>


namespace lf::dgfe {


using scalar_t = double;

class BoundingBox {
    public:
        BoundingBox(const lf::mesh::Entity &entity);
        /**
         * @brief maps from reference bounding box to the cell's bounding box
         * 
         * @param corners local points to be mapped into global coordinates
         * @return Eigen::Matrix Global points
         */
        Eigen::MatrixXd map(const Eigen::MatrixXd corners);

        /**
         * @brief Maps global coordinates into reference bounding box
         * 
         * @param corners global points
         * @return Eigen::MatrixXd local points
         */
        Eigen::MatrixXd inverseMap(const Eigen::MatrixXd corners);

        scalar_t det();

        /**
         * @brief returns entry (i,i) of the inverse jacobi
         */
        scalar_t inverseJacobi(unsigned i);

    private:
        Eigen::Vector2d translation_;
        Eigen::Matrix2d jacobi_;
        Eigen::Matrix2d inverse_jacobi_;
};

} //namespace lf::dgfe

#endif //LF_BOUNDING_BOX_H