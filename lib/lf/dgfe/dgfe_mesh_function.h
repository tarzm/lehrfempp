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
#include <lf/mesh/polytopic2d/polytopic2d.h>

#include "dgfe_space.h"
#include "dgfe_providers.h"
#include "integration.h"
#include "bounding_box.h"
#include "legendre_dgfe.h"

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
     * @brief evaluates the mesh function on a number of given global points. 
     * Global means that they will be mapped from the mesh into the reference bounding box inside the function.
     * 
     * @param e the polygon in which the points are located (inside or on the boundary)
     * @param global the global points on which the function is evaluated
     * @return std::vector<SCALAR> vector of function evaluations of the local points
     */
    std::vector<SCALAR> operator()(const lf::mesh::Entity& e, const Eigen::MatrixXd& global) const {
        auto dof_loc_coeffs = dgfe_space_->LocGlobMap().GlobalDofIndices(e);
        auto max_degree = dgfe_space_->MaxLegendreDegree();
        int n_points = global.cols();
        //initialize result vector with 0
        std::vector<SCALAR> result(n_points, 0.0);

        //get the bounding box of the polygon to map the points into the reference bounding box
        lf::dgfe::BoundingBox box(e);
        //result is calculated on the local coordinates
        auto local = box.inverseMap(global);

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

template<typename SCALAR, typename FUNCTOR> 
SCALAR L2ErrorSubTessellation(lf::dgfe::MeshFunctionDGFE<SCALAR> dgfe_function, FUNCTOR f, int max_degree){

    auto mesh_ptr = dgfe_function.Space()->Mesh();
    SCALAR error = 0.0;

    //lambda to calculate difference between exact solution and dgfe solution in one point
    auto errorAtPoint = [&dgfe_function, &f, max_degree](const lf::mesh::Entity *entity, Eigen::Vector2d global) -> SCALAR {
        LF_VERIFY_MSG(entity->RefEl() == lf::base::RefEl::kPolygon(), "Only implemented for polygons");
        
        lf::dgfe::BoundingBox box(*entity);
        return std::abs(dgfe_function(*entity, global)[0] - f(entity, global));
    };

    lf::dgfe::SubTessellationIntegrator<SCALAR, decltype(errorAtPoint)> integrator;

    //loop over cells
    for (auto cell : mesh_ptr->Entities(0)){
        error += integrator.integrate(cell, errorAtPoint, max_degree);
    }

    return error;
}


} //namespace lf::dgfe




#endif //MESH_FUNCTION_DGFE_H