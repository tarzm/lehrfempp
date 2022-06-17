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

DGFESpace::DGFESpace(std::shared_ptr<const lf::mesh::Mesh> mesh_p, size_type max_legendre_degree) : max_legendre_degree_(max_legendre_degree) ,
                     num_shape_funct_polygon_((max_legendre_degree_ == 1)? 4 : 9), dofh_(initializeDofHandler(mesh_p, num_shape_funct_polygon_)),
                     edge_polygon_adjacency_(initializeEdgePolygonAdjacency()) {
    LF_VERIFY_MSG(max_legendre_degree == 1 || max_legendre_degree == 2, "DGFE space currently only implemented for 1D legendre polynomials of degree 1 and 2 as basis functions");
}
                                        
const std::shared_ptr<const lf::mesh::Mesh> DGFESpace::Mesh() const {
    return dofh_.Mesh();
}

const lf::assemble::UniformDGFEDofHandler& DGFESpace::LocGlobMap() const {
    return dofh_;
}

size_type DGFESpace::NumRefShapeFunctions(const lf::mesh::Entity* entity){
    switch(entity->RefEl()){
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
    lf::base::RefEl ref_el = entity->RefEl();
    
    switch (ref_el){
        case lf::base::RefEl::kSegment():
            return edge_polygon_adjacency_(*entity);
        
        default:
            LF_VERIFY_MSG(false, "This functions only works for Segments");
    }
}

lf::assemble::UniformDGFEDofHandler DGFESpace::initializeDofHandler(std::shared_ptr<const lf::mesh::Mesh> mesh_p, size_type max_legendre_degree){
    //setup of dofhandler
    std::map<lf::base::RefEl, base::size_type> dofmap;
    dofmap.insert(std::pair<lf::base::RefEl, lf::base::size_type>(lf::base::RefEl::kPolygon(), num_shape_funct_polygon_));
    lf::assemble::UniformDGFEDofHandler dofhandler(mesh_p, dofmap);
    return dofhandler;
}

lf::mesh::utils::CodimMeshDataSet<lf::mesh::polytopic2d::PolygonPair> DGFESpace::initializeEdgePolygonAdjacency(){
    return lf::mesh::polytopic2d::EdgePolygonAdjacency(dofh_.Mesh());
}

size_type DGFESpace::MaxLegendreDegree() const{
    return max_legendre_degree_;
}



} //namespace lf::dgfe