/**
 * @file
 * @brief Class representing discrete space in the DGFE setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#ifndef DGFE_H
#define DGFE_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/assemble/assemble.h>



#include "legendre_dgfe.h"
#include "integration.h"

namespace lf::dgfe {

class DGFESpace {

public:

    DGFESpace(std::shared_ptr<const lf::mesh::Mesh> mesh_p, size_type max_legendre_degree);
    std::shared_ptr<const lf::mesh::Mesh> Mesh() const;
    const lf::assemble::DofHandler &LocGlobMap() const;
    size_type NumRefShapeFunctions(const lf::mesh::Entity*);
    lf::mesh::polytopic2d::PolygonPair AdjacentPolygons(const lf::mesh::Entity* entity);

private:
    /** maximum degree of 1D legendre polynomials used as basis fucntions*/
    size_type max_legendre_degree_;
    /** number of local shape functions on each polygon*/ 
    size_type num_shape_funct_polygon_;
    /** Local-to-global index map for the finite element space */
    lf::assemble::UniformDGFEDofHandler dofh_;
    /** data set which contains 2 (inner Segments) or 1 (boundary Segment) pointers to polygons
     * for each Segment*/
    lf::mesh::utils::CodimMeshDataSet<lf::mesh::polytopic2d::PolygonPair> edge_polygon_adjacency_;

};





} //namespace lf::dgfe

#endif //DGFE_O1_H