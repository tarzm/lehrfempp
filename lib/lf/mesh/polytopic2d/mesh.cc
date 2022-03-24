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

Mesh::Mesh(dim_t dim_world, NodeCoordList nodes, EdgeList edges, CellList cells, bool check_completeness)
           : dim_world_(dim_world) {

    //edges are not used at all! They are constructed from cells

    //if check_completenes == TRUE, we check for all Points to belong to a Segment. It is impossible
    //for Segments to not belong to a Polygon as they are only constructed in relation to a Polygon
    std::vector<bool> nodeHasSuperEntity;
    if (check_completeness) {
      nodeHasSuperEntity.resize(nodes.size(), false);
    }

    points_.reserve(nodes.size());
    polygons_.reserve(cells.size());

    //Global indexes of the entities
    size_type point_global_index = 0;
    size_type segment_global_index = 0;
    size_type polygon_global_index = 0;

    //map to check whether Segments already exist or not
    //maps the pair of Points of a Segment (first_node, second_node) to the index of the Segment
    std::map<const std::pair<size_type, size_type>, const size_type> SegmentIdxMap;
    using PointPairAndIdx = std::map<const std::pair<size_type, size_type>, const size_type>::value_type;
    
    //First fill the points_ vector of the mesh
    for (polytopic2d::Mesh::GeometryPtr &point_geo_ptr : nodes){
        points_.emplace_back(point_global_index, std::move(point_geo_ptr));
        point_global_index++;
    }

    //loop over edges to construct explicit edges
    for (std::pair<std::array<size_type, 2>, GeometryPtr> &edge : edges){
        //CHeck that edge does not already exist
        LF_VERIFY_MSG(SegmentIdxMap.find({edge.first[0], edge.first[1]}) == SegmentIdxMap.end()
                    && SegmentIdxMap.find({edge.first[1], edge.first[0]}) == SegmentIdxMap.end(),
                    "Somethign went wrong, an edge is passed to Mesh constructor that already exists.");

        //get pointers to endpoints
        const lf::mesh::hybrid2d::Point *p0_ptr = &points_[edge.first[0]];
        const lf::mesh::hybrid2d::Point *p1_ptr = &points_[edge.first[2]];

        //if geometry is nullptr construct it from endpoints
        GeometryPtr edge_geo_ptr(std::move(edge.second));
        if (edge_geo_ptr == nullptr){
            const Eigen::MatrixXd zero_point = base::RefEl::kPoint().NodeCoords();
            Eigen::Matrix<double, 2, 2> straight_edge_coords;
            straight_edge_coords.block<2, 1>(0, 0) = p0_ptr->Geometry()->Global(zero_point);
            straight_edge_coords.block<2, 1>(0, 1) = p1_ptr->Geometry()->Global(zero_point);
            edge_geo_ptr = std::make_unique<geometry::SegmentO1>(straight_edge_coords);
        }
        segments_.emplace_back(segment_global_index, std::move(edge_geo_ptr), p0_ptr, p1_ptr);

        segment_global_index++;
    }

    //First loop over cells to construct Segments
    for (auto cell_nodes_idx_vec : cells){
        int n_nodes = cell_nodes_idx_vec.size();

        //loop over nodes of the Polygon to construct the missing Segments
        for(int node_loc_idx = 0; node_loc_idx < n_nodes; node_loc_idx++){
            size_type current_node_glb_idx = cell_nodes_idx_vec.at(node_loc_idx);
            size_type next_node_glb_idx = cell_nodes_idx_vec.at((node_loc_idx + 1) % n_nodes);

            //check whether edge already exists, construct it if not
            if (SegmentIdxMap.find({current_node_glb_idx, next_node_glb_idx}) == SegmentIdxMap.end()
                    && SegmentIdxMap.find({next_node_glb_idx, current_node_glb_idx}) == SegmentIdxMap.end()) { //edge does not exist yet => construct it and add it to segments_

                const lf::mesh::hybrid2d::Point *current_node_ptr = &points_[current_node_glb_idx];
                const lf::mesh::hybrid2d::Point *next_node_ptr = &points_[next_node_glb_idx];

                //construct Geometry of Segment
                Eigen::Matrix<double, 2, 2> straight_edge_coords;
                const Eigen::MatrixXd zero_point = base::RefEl::kPoint().NodeCoords();
                straight_edge_coords.block<2, 1>(0, 0) = current_node_ptr->Geometry()->Global(zero_point);
                straight_edge_coords.block<2, 1>(0, 1) = next_node_ptr->Geometry()->Global(zero_point);
                GeometryPtr edge_geo_ptr = std::make_unique<geometry::SegmentO1>(straight_edge_coords);

                //Construct Segment and add it to segments_ vector
                segments_.emplace_back(segment_global_index, std::move(edge_geo_ptr), current_node_ptr, next_node_ptr);

                //Add it to the SegmentIdxMap
                SegmentIdxMap.insert(PointPairAndIdx({current_node_glb_idx, next_node_glb_idx}, segment_global_index));

                //increase global segment index
                segment_global_index++;
            }

            if(check_completeness){
              nodeHasSuperEntity[current_node_glb_idx] = true;
              nodeHasSuperEntity[next_node_glb_idx] = true;
            }

        } 
    }//now all Points and Segments of all the Polygons exist
    
    //Second loop over cells to construct the Polygons
    for (auto cell_nodes_idx_vec : cells){
        int n_nodes = cell_nodes_idx_vec.size();

        //vectors to be filled for the construction of the Polygon
        std::vector<const mesh::hybrid2d::Point*> corners(n_nodes, nullptr);
        std::vector<const mesh::hybrid2d::Segment*> edges(n_nodes, nullptr);

        size_type segment_loc_idx = 0;

        // loop over indexes of nodes
        for(int node_loc_idx = 0; node_loc_idx < n_nodes; node_loc_idx++){
            size_type current_node_glb_idx = cell_nodes_idx_vec.at(node_loc_idx);
            size_type next_node_glb_idx = cell_nodes_idx_vec.at((node_loc_idx + 1) % n_nodes);

            //add Point to corners vector
            const lf::mesh::hybrid2d::Point *current_node_ptr = &points_.at(current_node_glb_idx);
            LF_ASSERT_MSG(current_node_ptr->index() == current_node_glb_idx, "Node index does not coincide with its position in the points_ vector");
            corners.at(node_loc_idx) = current_node_ptr;

            //add Segment to edges vector
            if (SegmentIdxMap.find({current_node_glb_idx, next_node_glb_idx}) != SegmentIdxMap.end()){
                //edge exists in the order {current_node_glb_idx, next_node_glb_idx} and has positive local orientation
                auto edge_idx = SegmentIdxMap[{current_node_glb_idx, next_node_glb_idx}];
                edges.at(segment_loc_idx) = &segments_.at(edge_idx);
                segment_loc_idx++;
            } else if (SegmentIdxMap.find({next_node_glb_idx, current_node_glb_idx}) != SegmentIdxMap.end()){
                //edge exists in the order {next_node_glb_idx, current_node_glb_idx} and has negative local orientation
                auto edge_idx = SegmentIdxMap[{next_node_glb_idx, current_node_glb_idx}];
                edges.at(segment_loc_idx) = &segments_.at(edge_idx);
                segment_loc_idx++;
            } else {
                LF_ASSERT_MSG(false, "Something went wrong in the construction or indexing of the Segments");
            }
        }// END loop over indexes of nodes

        //contruct the Polygon
        polygons_.emplace_back(polygon_global_index, std::move(corners), std::move(edges));
        polygon_global_index++;
    
    } //end loop polygons 2

    //Set the entity_pointers_ vector's content
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

lf::mesh::utils::CodimMeshDataSet<PolygonPair> EdgePolygonAdjacency(std::shared_ptr<lf::mesh::Mesh> mesh_ptr){
    //initialize the CodimMeshDataSet
    PolygonPair init_pair = std::make_pair(nullptr, nullptr);
    lf::mesh::utils::CodimMeshDataSet<PolygonPair> adjacency(mesh_ptr, 1, init_pair);

    //set the values of the DataSet
    //loop over polygons
    for (const lf::mesh::Entity* polygon: mesh_ptr->Entities(0)){

        //loop over the polygons edges
        for(const lf::mesh::Entity* edge : polygon->SubEntities(1)){
            PolygonPair &polygon_pair = adjacency(*edge);
            if(adjacency(*edge).first == nullptr){
                adjacency(*edge).first = polygon;
            } else if (adjacency(*edge).second == nullptr){
                adjacency(*edge).second = polygon;
            } else {
                LF_VERIFY_MSG(false, "Something went wrong, both pointers of the edge polygon adjacency dataset are already set. An edge can only be adjacent to two polygons.\n");
            }
        }
    }  

    return adjacency;
}

} // namespace polytopic2d


