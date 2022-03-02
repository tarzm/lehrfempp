/**
 * @file
 * @brief Implementation of Polytopic2D Mesh class
 * @author Tarzis Maurer
 * @date   2022-02-26
 * @copyright ETH ZUrich
 */


#include "mesh.h"

#include <map>

namespace lf::mesh::polytopic2d {

nonstd::span<const Entity *const> Mesh::Entities(unsigned codim) const {
  LF_ASSERT_MSG(codim >= 0, "codim negative.");
  LF_ASSERT_MSG(codim <= dim_world_, "codim > dimWorld.");

  return entity_pointers_[codim];
}

Mesh::size_type Mesh::NumEntities(unsigned codim) const {
  switch (codim) {
    case 0:
      return polygons_.size();
    case 1:
      return segments_.size();
    case 2:
      return points_.size();
    default:
      LF_VERIFY_MSG(false, "codim out of bounds.");
  }
}

Mesh::size_type Mesh::NumEntities(lf::base::RefEl ref_el_type) const {
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return points_.size();
    }
    case lf::base::RefEl::kSegment(): {
      return segments_.size();
    }
    case lf::base::RefEl::kPolygon(): {
      return polygons_.size();
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal entity type");
      break;
    }
  }
  return 0;
}

Mesh::size_type Mesh::Index(const Entity &e) const {
  switch (e.Codim()) {
    case 0: {
      if (e.RefEl() == lf::base::RefEl::kPolygon()) {
        return dynamic_cast<const Polygon &>(e).index();
      }
      LF_VERIFY_MSG(false, "Illegal cell type");
    }
    case 1:
      return dynamic_cast<const lf::mesh::hybrid2d::Segment &>(e).index();
    case 2:
      return dynamic_cast<const lf::mesh::hybrid2d::Point &>(e).index();
    default:
      LF_VERIFY_MSG(false,
                    "Something is horribyl wrong, this entity has codim = " +
                        std::to_string(e.Codim()));
  }
}

const Entity *Mesh::EntityByIndex(dim_t codim, glb_idx_t index) const {
  LF_ASSERT_MSG(codim <= 2, "Illegal codimension " << codim);
  LF_ASSERT_MSG(index < NumEntities(codim),
                "Index " << index << " > " << NumEntities(codim));
  return entity_pointers_[codim][index];
}

bool Mesh::Contains(const Entity &e) const {
  switch (e.Codim()) {
    case 0:
      return !polygons_.empty() && &e >= &polygons_.front() && &e <= &polygons_.back();
    case 1:
      return &e >= &segments_.front() && &e <= &segments_.back();
    case 2:
      return &e >= &points_.front() && &e <= &points_.back();
    default:
      return false;
  }
}

Mesh::Mesh(dim_t dim_world, NodeCoordList nodes, EdgeList edges, CellList cells, bool check_completeness){
    
    //First fill the points_ vector of the mesh
    size_type node_index = 0;
    for (polytopic2d::Mesh::GeometryPtr &point_geo_ptr : nodes){
        points_.emplace_back(node_index, std::move(point_geo_ptr));
    }
}


} // namespace polytopic2d


