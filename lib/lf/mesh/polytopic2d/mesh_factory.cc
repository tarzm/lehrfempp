/**
 * @file
 * @brief Implementation of Polytopic2D MeshFactory class
 * @author Tarzis Maurer
 * @date   2022-03-01
 * @copyright ETH ZUrich
 */

#include "mesh_factory.h"


namespace lf::mesh::polytopic2d {


size_type MeshFactory::AddPoint(coord_t coord) {
  LF_ASSERT_MSG(coord.rows() == dim_world_,
                "coord has incompatible number of rows.");
  // Create default geometry object for a point from location vector
  polytopic2d::Mesh::GeometryPtr point_geo =
      std::make_unique<geometry::Point>(coord);
  nodes_.emplace_back(std::move(point_geo));
  return nodes_.size() - 1;
}

size_type MeshFactory::AddPoint(std::unique_ptr<geometry::Geometry>&& geometry){
    LF_ASSERT_MSG(false, "Add points with an Eigen::VectorXd as argument");
    return 0;
}

size_type MeshFactory::AddEntity(base::RefEl ref_el, const nonstd::span<const size_type>& nodes,
                                 std::unique_ptr<geometry::Geometry>&& geometry){
    LF_ASSERT_MSG(ref_el == lf::base::RefEl::kPolygon(), "RefEl must be of type kPolygon");
    std::vector<size_type> nodes_vec;
    for (const auto &node : nodes){
        LF_ASSERT_MSG(node < nodes_.size(), " Node " << node << " for " << ref_el.ToString() << "  must be inserted with AddPoint() first.");
        nodes_vec.emplace_back(node);
    }

    elements_.emplace_back(nodes_vec);
    return elements_.size() - 1;
}

std::shared_ptr<mesh::Mesh> MeshFactory::Build(){
    //Actual construction is done by Mesh object
    polytopic2d::Mesh::EdgeList edges_{};
    lf::mesh::polytopic2d::Mesh* mesh_ptr = new lf::mesh::polytopic2d::Mesh(dim_world_,
                         std::move(nodes_), std::move(edges_), std::move(elements_), check_completeness_);
    
    // Clear all information supplied to the MeshFactory object
    nodes_ = polytopic2d::Mesh::NodeCoordList{};  // .clear();
    edges_ = polytopic2d::Mesh::EdgeList{};       // .clear();
    elements_ = polytopic2d::Mesh::CellList{};    // clear();

    return std::shared_ptr<lf::mesh::polytopic2d::Mesh>(mesh_ptr);
}

} //namespace lf::mesh::polytopic2d