/**
 * @file
 * @brief discontinutity penalization functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
*/

#ifndef DISCONTINUITY_PENALIZATION_H
#define DISCONTINUITY_PENALIZATION_H

#include "integration.h"

namespace lf::dgfe {

using scalar_t = double;

/**
 * @brief returns a vector of areas which are the volumes of triangles made up of
 * the edge given as an argument and one other vertice. The ordering is counterclockwise,
 * such that the first value in the vector is the area of the triangle 
 * made up by endpoint_0, endpoint_1, and the next node of the polygon in counterclockwise direction
 * 
 * @param cell 
 * @return std::vector<scalar_t> 
 */
std::vector<scalar_t> simplexAreas(const lf::mesh::Entity &cell, size_type egde_loc_idx);








} //namespace lf::dgfe



#endif //DISCONTINUITY_PENALIZATION_H