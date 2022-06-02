/**
 * @file
 * @brief Implementation discrete space in the DGFE setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#include <lf/base/base.h>

#include "dgfe_space.h"

namespace lf::dgfe {

DGFESpace::DGFESpace(std::shared_ptr<const lf::mesh::Mesh> mesh_p, size_type max_legendre_degree){
    LF_VERIFY_MSG(max_legendre_degree == 1 || max_legendre_degree == 2, "DGFE space currently only implemented for 1D legendre polynomials of degree 1 and 2");
    max_legendre_degree_ = max_legendre_degree;

    if (max_legendre_degree_ == 1){
        num_shape_funct_polygon_ = 4;
    } else {
        num_shape_funct_polygon_ = 9;
    }

    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), num_shape_funct_polygon_));
    dofh_ = lf::assemble::UniformDGFEDofHandler(std::move(mesh_p), dofmap);

    //setup of adjacency data set
    edge_polygon_adjacency_ = lf::mesh::polytopic2d::EdgePolygonAdjacency(dofh_.Mesh());
}
                                        
std::shared_ptr<const lf::mesh::Mesh> DGFESpace::Mesh() const {
    return dofh_.Mesh();
}

const lf::assemble::DofHandler& DGFESpace::LocGlobMap() const {
    return dofh_;
}

size_type DGFESpace::NumRefShapeFunctions(const lf::mesh::Entity* entity){
    switch(entity.RefEl()){
        case lf::base::RefEl::kPolygon():
            return num_shape_funct_polygon_;

        case lf::base::RefEl::kSegment():
        case lf::base::RefEl::kPoint():
            return 0;
        
        default:
            LF_ASSERT_MSG(false, "Invalid type of entity. DGFE space only accepts Points, Segments and Polygons");

    }
}

lf::mesh::polytopic2d::PolygonPair DGFESpace::AdjacentPolygons(const lf::mesh::Entity* entity){
    lf::base::RefEl ref_el = entity.RefEL();
    
    switch (ref_el){
        case lf::base::RefEl::kSegment():
            return edge_polygon_adjacency_(entity);
        
        default:
            LF_VERIFY_MSG(false, "This functions only works for Segments");
    }
}





} //namespace lf::dgfe