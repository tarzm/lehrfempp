#include "mesh_factory.h"
#include "entity.h"

namespace lf::mesh::hybrid2d {

MeshFactory::size_type MeshFactory::AddPoint(coord_t coord) {
  LF_ASSERT_MSG(!built_, "Build() already called.");
  LF_ASSERT_MSG(coord.rows() == dim_world_,
                "coord has incompatible number of rows.");
  nodes_.emplace_back(std::move(coord));
  return nodes_.size() - 1;
}

MeshFactory::size_type MeshFactory::AddPoint(
    std::unique_ptr<geometry::Geometry>&& /*geometry*/) {
  LF_ASSERT_MSG(false, "Not implemented");
  return -1;
}

MeshFactory::size_type MeshFactory::AddEntity(
    base::RefEl ref_el, const base::ForwardRange<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) {
  LF_ASSERT_MSG(!built_, "Build() already called.");
  LF_ASSERT_MSG(ref_el.Dimension() > 0,
                "Use AddNode() to add a node to a mesh.");
  LF_ASSERT_MSG(ref_el.Dimension() <= 2, "ref_el.Dimension > 2");
  if (geometry != nullptr) {
    LF_ASSERT_MSG(geometry->DimGlobal() == dim_world_,
                  "geometry->DimGlobal() != dim_world_");
    LF_ASSERT_MSG(geometry->RefEl() == ref_el, "ref_el != geometry->RefEl()");
  }

  if (ref_el == base::RefEl::kSegment()) {
    std::array<size_type, 2> ns{{0, 0}};
    unsigned char count = 0;
    for (size_type n : nodes) {
      LF_ASSERT_MSG(n < nodes_.size(),
                    "node " << n
                            << " specified in call to AddEntity must be "
                               "inserted with AddNode() first.");
      LF_ASSERT_MSG(
          count < 2,
          "ref_el = segment, but nodes contains more than 2 node indices");
      ns[count] = n;
      ++count;
    }
    LF_ASSERT_MSG(count == 2,
                  "ref_el = segment but size of nodes was " << count);
    edges_.emplace_back(ns, std::move(geometry));
    return edges_.size() - 1;
  }

  // otherwise its an element:
  LF_ASSERT_MSG(geometry, "Geometry is required for elements (codim=0)");
  std::array<size_type, 4> ns{};
  unsigned char count = 0;
  for (auto& n : nodes) {
    LF_ASSERT_MSG(n < nodes_.size(),
                  "node " << n
                          << " specified in call to AddEntity must be "
                             "inserted with AddNode() first.");
    LF_ASSERT_MSG(count < ref_el.NumNodes(),
                  "ref_el = " << ref_el << ", but nodes contains " << count + 1
                              << "node indices");
    ns[count] = n;
    ++count;
  }
  LF_ASSERT_MSG(count == ref_el.NumNodes(),
                "ref_el.NumNodes() = " << ref_el.NumNodes()
                                       << ", but argument nodes contained "
                                       << count << " nodes");

  // if we have added a triangle, set the fourth node to -1.
  if (count == 3) {
    ns[3] = -1;
  }
  elements_.emplace_back(ns, std::move(geometry));
  return elements_.size() - 1;
}

std::shared_ptr<mesh::Mesh> MeshFactory::Build() {
  built_ = true;
  return std::shared_ptr<Mesh>(new Mesh(
      dim_world_, std::move(nodes_), std::move(edges_), std::move(elements_)));
}

}  // namespace lf::mesh::hybrid2d