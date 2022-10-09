/**
 * @file
 * @brief Implementation of Polytopic2D MeshFactory class
 * @author Tarzis Maurer
 * @date   2022-03-01
 * @copyright ETH ZUrich
 */

#include <cmath>

#include "mesh_factory.h"

#define COORD_TOLERANCE 1e-7

namespace lf::mesh::polytopic2d {


size_type MeshFactory::AddPoint(coord_t coord) {
    LF_ASSERT_MSG(coord.rows() == dim_world_,
                "coord has incompatible number of rows.");

    double x = coord[0];
    double y = coord[1];
    if (unit_square_){
        double x = coord[0];
        double y = coord[1];
        if (std::abs(x) < COORD_TOLERANCE){
            coord[0] = 0.0;
        }
        if (std::abs(x - 1.0) < COORD_TOLERANCE){
            coord[0] = 1.0;
        }
        if (std::abs(y) < COORD_TOLERANCE){
            coord[1] = 0.0;
        }
        if (std::abs(y - 1.0) < COORD_TOLERANCE){
            coord[1] = 1.0;
        }
    }
    
    // Create default geometry object for a point from location vector
    polytopic2d::Mesh::GeometryPtr point_geo = std::make_unique<geometry::Point>(coord);
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
    mesh::Mesh* mesh_ptr = new lf::mesh::polytopic2d::Mesh(dim_world_,
                         std::move(nodes_), std::move(elements_), check_completeness_);
    
    // Clear all information supplied to the MeshFactory object
    nodes_ = polytopic2d::Mesh::NodeCoordList{};  // .clear();
    elements_ = polytopic2d::Mesh::CellList{};    // clear();

    return std::shared_ptr<mesh::Mesh>(mesh_ptr);
}


std::shared_ptr<lf::mesh::Mesh> polytopicFromHybrid2D(std::shared_ptr<const lf::mesh::Mesh> mesh_ptr){
    using coord_t = Eigen::Vector2d;
    using size_type = lf::mesh::Mesh::size_type;

    // Obtain mesh factory
    std::unique_ptr<lf::mesh::polytopic2d::MeshFactory> mesh_factory_ptr = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2);

    //Add all nodes
    for (auto point : mesh_ptr->Entities(2)){
        auto coord = lf::geometry::Corners(*(point->Geometry()));
        mesh_factory_ptr->AddPoint(coord);
    }

    //Add all polygons
    for (auto cell : mesh_ptr->Entities(0)){
        std::vector<unsigned> node_indices;
        //loop over nodes of cell
        for (auto point : cell->SubEntities(2)){
            node_indices.push_back(mesh_ptr->Index(*point));
        }
        if(cell->RefEl() == lf::base::RefEl::kTria()){
            mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type, 3>{node_indices[0], node_indices[1], node_indices[2]}, nullptr);
        } else { //Quad
            mesh_factory_ptr->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type, 4>{node_indices[0], node_indices[1], node_indices[2], node_indices[3]}, nullptr);
        }
        
    }

    return mesh_factory_ptr->Build();
}

} //namespace lf::mesh::polytopic2d