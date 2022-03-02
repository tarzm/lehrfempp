/** @file
 *  @brief Entity implementation for polygons
 *  @author Tarzis Maurer
 *  @date 2022-02-14
 *  @copyright ETH Zurich
 */

#ifndef POLYGON_H
#define POLYGON_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

namespace lf::mesh::polytopic2d{

/**
 * @brief Describes a polygonal cell for a 2D polytopic mesh
 *
 * A polytopic cell is defined by ordered vectors of references to its nodes
 * and its edges; internal consistency is required
 * @note Every `Segment` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Polygon : public mesh::Entity {
    using size_type = mesh::Mesh::size_type;

    public:
        /**
         * @brief constructor, is called from MeshFactory
         * @param index index of the entity to be created; will usually be
         * retrieved via the `Index()` method of `Mesh`
         * @param corners vector of pointers to nodes of type mesh::hybrid2d::Point
         * @param edges vector of pointers to edges of type mesh::hybrid::Segment
         * 
         * The ordering of the sub-entities in the vectors is such that edge i connects node i and i+1
         *
         */
        explicit Polygon(size_type index,
                        std::vector<const mesh::hybrid2d::Point*> corners,
                        std::vector<const mesh::hybrid2d::Segment*> edges);
    
        /** @brief an edge is an entity of co-dimension 1 */
        [[nodiscard]] unsigned Codim() const override { return 0; }

        /**
         * @brief returns the number of nodes of the polygon
         */
        size_type NumNodes();

        /** @brief access to index of an entity */
        [[nodiscard]] size_type index() const { return index_; }

        /**
         * @brief returns the nodes of the polygion in matrix format
         */
        Eigen::MatrixXd Corners();



    private:
        size_type index_ = -1;                                      // zero-based index of this entity.
        std::unique_ptr<geometry::Geometry> geometry_ = nullptr;    // shape information => not used in polytopic setting
        std::vector<const mesh::hybrid2d::Point*> nodes_{};         // nodes = corners of cell
        std::vector<const mesh::hybrid2d::Segment*> edges_{};       // edges of the cells
        std::vector<lf::mesh::Orientation> edge_ori_{};             // orientation of edges (set in constructor)
        Entity* this_ = nullptr;                                    // needed for SubEntity()
        Eigen::MatrixXd corners_;
};

} //namespace lf::mesh::hybrid2d

#endif // POLYGON_H