/**
 * @file
 * @brief Functionalities for decomposing a polygon into geometry objects of triangles for integration
 * @author Tarzis Maurer
 * @date June
 * @copyright ETH Zurich
*/

#ifndef SUB_TESSELLATION_H
#define SUB_TESSELLATION_H

#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/geometry/geometry.h>


namespace lf::dgfe {

std::vector<std::unique_ptr<lf::geometry::TriaO1>> subTessellation(const lf::mesh::Entity *polygon);


} //namespace lf::dgfe



#endif // SUB_TESSELLATION_H