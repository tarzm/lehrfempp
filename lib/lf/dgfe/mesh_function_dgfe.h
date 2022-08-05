/**
 * @file
 * @brief Class representing mesh functions in the DGFE setting
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef MESH_FUNCTION_DGFE_H
#define MESH_FUNCTION_DGFE_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>

#include "dgfe_space.h"
#include "dgfe_providers.h"
#include "integration.h"
#include "bounding_box.h"
#include "legendre_dgfe.h"
#include "mesh_function_global.h"

namespace lf::dgfe {


template<typename SCALAR>
class MeshFunctionDGFE {

public:

    /**
     * @brief Create a new mesh function in the DGFE setting
     * 
     * @param dgfe_space approximation space in which the function lies
     * @param coeff_vector coefficients of the basis expansion of the DGFE space
     */
    MeshFunctionDGFE(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> coeff_vector)
                    : dgfe_space_(std::move(dgfe_space)), dof_vector_(std::move(coeff_vector)), num_shape_funct_polygon_((dgfe_space_->MaxLegendreDegree() == 1)? 4 : 9){}

    /**
     * @brief evaluates the mesh function on a number of given local points. 
     * local means that they have been mapped from the bounding box of the polygon into the referenc ebounding box before
     * 
     * @param e the polygon in which the points are located (inside or on the boundary)
     * @param local the local points on which the function is evaluated
     * @return std::vector<SCALAR> vector of function evaluations of the local points
     */
    std::vector<SCALAR> operator()(const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const {
        auto dof_loc_coeffs = dgfe_space_->LocGlobMap().GlobalDofIndices(e);
        auto max_degree = dgfe_space_->MaxLegendreDegree();
        int n_points = local.cols();
        //initialize result vector with 0
        std::vector<SCALAR> result(n_points, 0.0);

        //loop over basis functions and dof coefficients
        for (int i = 0; i < num_shape_funct_polygon_; i++){
            auto degrees = lf::dgfe::multiIndexToDegree(i, max_degree);
            auto legendre_degree_x = degrees.first;
            auto legendre_degree_y = degrees.second;

            //loop over local points
            for (int loc_idx = 0; loc_idx < n_points; loc_idx++){
                result[loc_idx] += dof_vector_[dof_loc_coeffs[i]] * legendre_polynomial(legendre_degree_x, local(0, loc_idx)) * legendre_polynomial(legendre_degree_y, local(1, loc_idx));
            }
        }
        return result;
    }

    /**
     * @brief getter method for the DGFE space
     */
    std::shared_ptr<const lf::dgfe::DGFESpace> Space(){
        return dgfe_space_;
    }


private:

    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> dof_vector_;
    /** number of local shape functions on each polygon*/ 
    size_type num_shape_funct_polygon_; 
};


template<typename SCALAR>
class MeshFunctionGradDGFE {

public:

    /**
     * @brief Create a new gradient mesh function in the DGFE setting
     * 
     * @param dgfe_space approximation space in which the function lies
     * @param coeff_vector coefficients of the basis expansion of the DGFE space
     */
    MeshFunctionGradDGFE(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> coeff_vector)
                    : dgfe_space_(std::move(dgfe_space)), dof_vector_(std::move(coeff_vector)), num_shape_funct_polygon_((dgfe_space_->MaxLegendreDegree() == 1)? 4 : 9){}

    /**
     * @brief evaluates the mesh function on a number of given local points. 
     * local means that they have been mapped from the bounding box of the polygon into the reference bounding box before
     * 
     * @param e the polygon in which the points are located (inside or on the boundary)
     * @param local the local points on which the function is evaluated
     * @return std::vector<Eigen::Matrix<SCALAR, 2, 1>> vector of function evaluations of the local points
     */
    std::vector<Eigen::Matrix<SCALAR, 2, 1>> operator()(const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const {
        auto dof_loc_coeffs = dgfe_space_->LocGlobMap().GlobalDofIndices(e);
        auto max_degree = dgfe_space_->MaxLegendreDegree();
        int n_points = local.cols();
        //initialize result vector with 0
        std::vector<Eigen::Matrix<SCALAR, 2, 1>> result(n_points, Eigen::Matrix<SCALAR, 2, 1>::Zero());

        lf::dgfe::BoundingBox box(e);
        
        //loop over basis functions and dof coefficients
        //starts at 1 because derivative of a constant function is always 0
        for (int basis = 1; basis < num_shape_funct_polygon_; basis++){

            //loop over local points
            for (int loc_idx = 0; loc_idx < n_points; loc_idx++){
                Eigen::Matrix<SCALAR, 2, 1> nabla_basis{    legendre_basis_dx(basis, dgfe_space_->MaxLegendreDegree(), local.col(loc_idx)) * box.inverseJacobi(0),
                                                            legendre_basis_dy(basis, dgfe_space_->MaxLegendreDegree(), local.col(loc_idx)) * box.inverseJacobi(1)};

                result[loc_idx] += dof_vector_[dof_loc_coeffs[basis]] * nabla_basis;
            }
        }
        return result;
    }

    /**
     * @brief getter method for the DGFE space
     */
    std::shared_ptr<const lf::dgfe::DGFESpace> Space(){
        return dgfe_space_;
    }


private:

    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> dof_vector_;
    /** number of local shape functions on each polygon*/
    size_type num_shape_funct_polygon_; 
};

} //namespace lf::dgfe

#endif //MESH_FUNCTION_DGFE_H