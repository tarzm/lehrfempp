/** @file
 *  @brief Implementation of Polygon class
 *  @author Tarzis Maurer
 *  @date 2022-02-14
 *  @copyright ETH Zurich
 */

#include "polygon.h"

namespace lf::mesh::polytopic2d {

    Polygon::Polygon(size_type index,
                    std::vector<const mesh::hybrid2d::Point*> nodes,
                    std::vector<const mesh::hybrid2d::Segment*> edges)
        : index_(index),
        nodes_(nodes),
        edges_(edges),
        edge_ori_(),
        this_(this){
        
        LF_VERIFY_MSG(nodes.size() == edges.size(), "Number of nodes and edges are not equal");

        for (const auto node : nodes){
            LF_VERIFY_MSG(node != nullptr, "Invalid pointer to a corner");
        }

        for (const auto edge : edges){
            LF_VERIFY_MSG(edge != nullptr, "Invalid pointer to an edge");
        }
        
        size_type num_nodes = nodes.size();

        //constrution of corners_
        corners_.resize((nodes.at(0)->RefEl()).Dimension(), num_nodes);
        for (int node_loc_idx = 0; node_loc_idx < num_nodes; node_loc_idx++){
            corners_.col(node_loc_idx) = lf::geometry::Corners(*(nodes.at(node_loc_idx)->Geometry()));
        }
        
        // Finally set relative orientations for the edges. Edge i has positive
        // orientation, if its first node agrees with vertex i
        for (int ed_loc_idx = 0; ed_loc_idx < num_nodes; ed_loc_idx++) {
            // Fetch nodes of current edge
            auto ed_nodes = edges_.at(ed_loc_idx)->SubEntities(1);
            edge_ori_.at(ed_loc_idx) = (ed_nodes[0] == nodes_.at(ed_loc_idx))
                                        ? lf::mesh::Orientation::positive
                                        : lf::mesh::Orientation::negative;
        }
    }

    mesh::Mesh::size_type Polygon::NumNodes(){
        return nodes_.size();
    }

    Eigen::MatrixXd Polygon::Corners(){
        return corners_;
    }

    nonstd::span<const Entity* const> Polygon::SubEntities(unsigned rel_codim) const {
        auto l = [&](auto i) -> const mesh::Entity& { return **i; };
        switch (rel_codim){
            case 2:
                return {reinterpret_cast<const Entity* const*>(&corners_[0]), NumNodes();}
            case 1:
                return {reinterpret_cast<const Entity* const*>(&edges_[0]), NumNodes();}
            case 0:
                return {&this_, 1};
            default:
                LF_VERIFY_MSG(false, "Triangle: rel_codim out of range");
        }
    }





}   //  namespace lf::mesh::polytopic2d