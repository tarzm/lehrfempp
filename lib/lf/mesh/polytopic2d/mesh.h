/**
 * @file
 * @brief Declares the Polytopic2d Mesh class
 * @author Tarzis Maurer
 * @date   2022-02-26
 * @copyright ETH Zurich
 */

#ifndef MESH_POLYTOPIC2D_H
#define MESH_POLYTOPIC2D_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <Eigen/Eigen>

#include "polygon.h"

namespace lf::mesh::polytopic2d {

using size_type = lf::base::size_type;
using dim_t = lf::base::dim_t;
using sub_idx_t = lf::base::sub_idx_t;
using glb_idx_t = lf::base::glb_idx_t;

class MeshFactory;

/** @brief Polytopic 2D mesh type compliant with abstract mesh interface
 * 
 */
class Mesh: public mesh::Mesh {
    public:
        [[nodiscard]] unsigned DimMesh() const override { return 2; }
        [[nodiscard]] unsigned DimWorld() const override { return dim_world_; }
        [[nodiscard]] nonstd::span<const Entity* const> Entities(unsigned codim) const override;
        [[nodiscard]] size_type NumEntities(unsigned codim) const override;
        [[nodiscard]] size_type NumEntities(lf::base::RefEl ref_el_type) const override;
        [[nodiscard]] size_type Index(const Entity& e) const override;
        [[nodiscard]] const mesh::Entity* EntityByIndex(dim_t codim, glb_idx_t index) const override;
        [[nodiscard]] bool Contains(const mesh::Entity& e) const override;
    

    private:
        dim_t dim_world_{};
        /** @brief array of 0-dimensional entity object of co-dimension 2 */
        std::vector<hybrid2d::Point> points_;
        /** @brief array of 1-dimensional entity object of co-dimension 1 */
        std::vector<hybrid2d::Segment> segments_;
        /** @brief array of polygonal cell objects, co-dimension 0 */
        std::vector<polytopic2d::Polygon> polygons_;
    

        /** @brief Auxiliary array of cell (co-dim ==0 entities) pointers
         *
         * This array serves two purposes. It facilitates the construction
         * of a range covering all the cells. It is also required to retrieving
         * the entity of a specific codimension belonging to an index.
         */
        std::array<std::vector<const mesh::Entity*>, 3> entity_pointers_;

        /** @brief Data types for passing information about mesh entities */
        using GeometryPtr = std::unique_ptr<geometry::Geometry>;
        using NodeCoordList = std::vector<GeometryPtr>;
        using EdgeList = std::vector<std::pair<std::array<size_type, 2>, GeometryPtr>>;
        using CellList = std::vector<std::vector<size_type>>;

        /**
         * @brief Construction of mesh from information gathered in a MeshFactory
         * @param dim_world Dimension of the ambient space.
         * @param nodes sequential container of node coordinates
         * @param edges sequential container of pairs of
         *               (i) vectors of indices of the nodes of an edge
         *               (ii) pointers to the geometry object describing an edge
         * @param cells sequential container of pairs of
         *               (i) vectors of indices of the nodes of a cell
         *               (ii) pointers to the geometry object for the cell. 
         *                    Here in the polytopic case those are nullptrs.
         * @param check_completeness If set to true, the constructor will check that
         * the mesh is topologically complete. That means that every entity with
         * codimension `codim>0` is a subentity of at least one entity with
         * codimension `codim-1`. If `check_completeness = true` and the mesh is not
         * complete, an assert will fail.
         *
         * ### Shape guessing
         *
         * An extreme situation is marked by passing node positions as the only
         * geometric information together with topological node-cell incidence
         * relationships. In this case the constructor will build a mesh with straight
         * edges throughout, type lf::geometry::SegmentO1.
         *
         * ### Missing entities
         *
         * @note all edges are currently inserted based on the information about the cells.
         *
         * @note the position of node information the `nodes` array and of cell
         *        information in the `cells` array, respectively,
         *        determines the interpretation of the index numbers,
         *        that is the n-th node in the container has index n-1.
         *
         */
        Mesh(dim_t dim_world, NodeCoordList nodes, CellList cells, bool check_completeness);

        friend class MeshFactory;
};

using PolygonPair = std::pair<const lf::mesh::Entity*, const lf::mesh::Entity*>;


/**
 * @brief Constructs A CodimMeshDataSet that contains the adjacencies of the Segements. Each Segment is adjacent to either two Polygons
 *        (inner Segment) or one Polygon (boundary Segment).
 * 
 * @param mesh_ptr The mesh used.
 * @return lf::mesh::utils::CodimMeshDataSet<PolygonPair> The constructed CodimMeshDataSet
 */
lf::mesh::utils::CodimMeshDataSet<PolygonPair> EdgePolygonAdjacency(std::shared_ptr<lf::mesh::Mesh> mesh_ptr);

}  // namespace lf::mesh::hybrid2d

#endif // MESH_POLYTOPIC2D_H