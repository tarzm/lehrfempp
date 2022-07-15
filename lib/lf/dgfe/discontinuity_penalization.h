/**
 * @file
 * @brief discontinutity penalization functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
*/

#ifndef DISCONTINUITY_PENALIZATION_H
#define DISCONTINUITY_PENALIZATION_H

#include "dgfe_space.h"
#include "integration.h"

namespace lf::dgfe {

using scalar_t = double;

/**
 * @brief returns a vector of areas which are the volumes of triangles made up of
 * the edge given as an argument and one other vertice. The ordering is counterclockwise,
 * such that the first value in the vector is the area of the triangle 
 * made up by endpoint_0, endpoint_1, and the next node of the polygon in counterclockwise direction
 * 
 * @note this only works for convex polygons!
 */
std::vector<scalar_t> simplexAreas(const lf::mesh::Entity &cell, size_type egde_loc_idx);



class DiscontinuityPenalization{
    
    public:
        DiscontinuityPenalization(std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr, scalar_t c_inv_constant, scalar_t c_sigma_constant) : 
                                dgfe_space_ptr_(std::move(dgfe_space_ptr)), c_inv_const_(c_inv_constant), c_sigma_const_(c_sigma_constant) {}

        scalar_t operator()(const lf::mesh::Entity &edge, scalar_t A_f) const ;

        std::shared_ptr<const lf::dgfe::DGFESpace> dgfeSpace();
    

    private:
        std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_ptr_;
        scalar_t c_inv_const_;
        scalar_t c_sigma_const_;


};



} //namespace lf::dgfe



#endif //DISCONTINUITY_PENALIZATION_H