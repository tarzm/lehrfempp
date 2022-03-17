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

enum class EdgeExistence : int {positive = 1, negative = -1, inexistent = 0};

Mesh::Mesh(dim_t dim_world, NodeCoordList nodes, CellList cells, bool check_completeness) : dim_world_(dim_world){

    //edges are not used at all! They are constructed from cells

    //if check_completenes == TRUE, we check for all Points to belong to a Segment. It is impossible
    //for Segments to not belong to a Polygon as they are only constructed in relation to a Polygon
    std::vector<bool> nodeHasSuperEntity;
    if (check_completeness) {
      nodeHasSuperEntity.resize(nodes.size(), false);
    }

    points_.reserve(nodes.size());
    polygons_.reserve(cells.size());
    
    //First fill the points_ vector of the mesh
    size_type node_index = 0;
    for (polytopic2d::Mesh::GeometryPtr &point_geo_ptr : nodes){
        points_.emplace_back(node_index, std::move(point_geo_ptr));
        //DEGBUGGING:
        node_index++;
    }

    //|||||||||||||||||||||||||Construct all Segments from the cells vector then construct Polygons
    size_type segment_index = 0;
    size_type polygon_index = 0;

    //map to check whether edges already exist or not
    //maps the pair of nodes of an edge (first_node, second_node) to the index of the edge
    std::map<const std::pair<size_type, size_type>, const size_type> EdgeIdxMap;
    using EdgePairAndIdx = std::map<const std::pair<size_type, size_type>, const size_type>::value_type;

    //loop over Polygons
    for (auto cell_nodes_idx : cells){
        int n_nodes = cell_nodes_idx.size();

        //vectors to be filled for construction of each Polygon
        std::vector<const mesh::hybrid2d::Point*> corners;
        std::vector<size_type> edge_indexes;

        size_type edges_current_size = segments_.size();
        segments_.reserve(edges_current_size + n_nodes);


        //loop over nodes of the Polygon to construct Edges first, then Polygons
        for(int node_idx = 0; node_idx < n_nodes; node_idx++){
            size_type current_node = cell_nodes_idx.at(node_idx);
            size_type next_node = cell_nodes_idx.at((node_idx + 1) % n_nodes);

            //add Point to corners vector
            const lf::mesh::hybrid2d::Point *node_ptr = &points_.at(current_node);
            LF_ASSERT_MSG(node_ptr->index() == current_node, "Node index does not coincide with its position in the points_ vector");
            corners.push_back(node_ptr);

            //index of the next edge that will later be added to edge_indexes;
            size_type next_segment_idx;

            //check whether edge already exists, construct it if necessary
            //and set next_segment_idx correctly
            if (EdgeIdxMap.find({current_node, next_node}) != EdgeIdxMap.end()){
                //edge exists in the order {current_node, next_node} => add it to edges
                auto edge_idx = EdgeIdxMap[{current_node, next_node}];
                next_segment_idx = edge_idx;
            } else if (EdgeIdxMap.find({next_node, current_node}) != EdgeIdxMap.end()){
                //edge exists in the order {next_node, current_node} => add it to edges
                auto edge_idx = EdgeIdxMap[{next_node, current_node}];
                next_segment_idx = edge_idx;
            } else { //edge does not exist yet => construct it, add it to segments_ and edges
                //retrieve pointer to next node (current node is defined above as node_ptr)
                const lf::mesh::hybrid2d::Point *next_node_ptr = &points_[next_node];
                //construct Geometry of Segment
                Eigen::Matrix<double, 2, 2> straight_edge_coords;
                const Eigen::MatrixXd zero_point = base::RefEl::kPoint().NodeCoords();
                straight_edge_coords.block<2, 1>(0, 0) = node_ptr->Geometry()->Global(zero_point);
                straight_edge_coords.block<2, 1>(0, 1) = next_node_ptr->Geometry()->Global(zero_point);
                GeometryPtr edge_geo_ptr = std::make_unique<geometry::SegmentO1>(straight_edge_coords);
                //Construct Segment and add it to segments_ vector
                segments_.emplace_back(segment_index, std::move(edge_geo_ptr), node_ptr, next_node_ptr);
                //Add it to the EdgeIdxMap
                EdgeIdxMap.insert(EdgePairAndIdx({current_node, next_node}, segment_index));
                //Set next_segment_idx to current one
                next_segment_idx = segment_index;
                segment_index++;

            }
            if(check_completeness){
              nodeHasSuperEntity[current_node] = true;
              nodeHasSuperEntity[next_node] = true;
            }

            //Add next_segment_idx to edge_indexes
            edge_indexes.push_back(next_segment_idx);

        } //now all nodes and edges of the Polygon exist => construct the Polygon

        //vector of Segments that will be passed to constructor of Polygon
        std::vector<const mesh::hybrid2d::Segment*> edges;
        for (int i = 0; i < edge_indexes.size(); i++){
            const lf::mesh::hybrid2d::Segment *edge_ptr = &(segments_[edge_indexes[i]]);
            edges.push_back(std::move(edge_ptr));
        }

        //construct Polygon
        polygons_.emplace_back(polygon_index, std::move(corners), std::move(edges));

        polygon_index++;
    }
    
    //Fill the entity_pointers_ array
    entity_pointers_[0] = std::vector<const lf::mesh::Entity *>(polygons_.size(), nullptr);
    entity_pointers_[1] = std::vector<const lf::mesh::Entity *>(segments_.size(), nullptr);
    entity_pointers_[2] = std::vector<const lf::mesh::Entity *>(points_.size(), nullptr);
    
    
    for (auto &pol : polygons_) {
        entity_pointers_[0].at(pol.index()) = &pol;
    }
    for (auto &s : segments_) {
        entity_pointers_[1].at(s.index()) = &s;
    }
    for (auto &p : points_) {
        entity_pointers_[2].at(p.index()) = &p;
    }
    
    
    if (check_completeness) {
        // Check that all nodes have a super entity:
        for (std::size_t i = 0; i < nodes.size(); ++i) {
        LF_VERIFY_MSG(nodeHasSuperEntity[i],
                        "Mesh is incomplete: Node with global index "
                            << i << " is not part of any edge.");
        }
    }
}

} // namespace polytopic2d


