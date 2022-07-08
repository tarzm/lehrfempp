/**
 * @file
 * @brief Functionalities inlcuding average and jump operators
 * which are present in the general DGFEM form of linear degenerate
 * second order convection-diffusion-reaction bundary value problem
 * 
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef DGFE_AUX_OPERATORS_H
#define DGFE_AUX_OPERATORS_H


#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>


namespace lf::dgfe {


template <typename SCALAR, typename MESHFUNCTION>
SCALAR averageOperator(const lf::dgfe::DGFESpace *dgfe_space_ptr, const lf::mesh::Entity &entity, MESHFUNCTION f, const Eigen::Vector2d &local){
    LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kSegment(), "Only implemented for Segments");
    //get pointers to adjacent polygons
    auto polygonPair = dgfe_space_ptr->AdjacentPolygons(entity);
    auto polygon_1 = polygonPair.first;
    auto polygon_2 = polygonPair.second;
    return 0.5 * (f(polygon_1, local) + f(polygon_2, local));
}


} //namespace lf::dgfe


#endif //DGFE_AUX_OPERATORS_H