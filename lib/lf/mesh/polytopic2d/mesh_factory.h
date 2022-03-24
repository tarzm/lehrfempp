/**
 * @file
 * @brief Declares the Polytopic2d MeshFactory class
 * @author Tarzis Maurer
 * @date   2022-03-01
 * @copyright ETH Zurich
 */

#ifndef MESH_FACTORY_POLYTOPIC2D_H
#define MESH_FACTORY_POLYTOPIC2D_H

#include <lf/mesh/mesh.h>

#include "mesh.h"

namespace lf::mesh::polytopic2d {

class MeshFactory : public mesh::MeshFactory{
    public:
        MeshFactory(const MeshFactory&) = delete;
        MeshFactory(MeshFactory&&) = delete;
        MeshFactory& operator=(const MeshFactory&) = delete;
        MeshFactory& operator=(MeshFactory&&) = delete;

        /**
         * @brief Construct a new builder that can be used to construct a new polytopic2d
         *        mesh.
         * @param dim_world The dimension of the euclidean space in which the
         *                  mesh is embedded.
         * @param check_completeness If set to true, calling Build() will check that
         * the mesh is topologically complete. That means that every entity with
         * codimension `codim>0` is a subentity of at least one entity with
         * codimension `codim-1`. If `check_completeness = true` and the mesh is not
         * complete, an assert will fail.
         */
        explicit MeshFactory(dim_t dim_world, bool check_completeness = true)
            : dim_world_(dim_world), check_completeness_(check_completeness) {}

        
        [[nodiscard]] dim_t DimWorld() const override { return dim_world_; }

        [[nodiscard]] dim_t DimMesh() const override { return 2; }

        // NOLINTNEXTLINE(modernize-use-nodiscard)
        size_type AddPoint(coord_t coord) override;

        // NOLINTNEXTLINE(modernize-use-nodiscard)
        size_type AddPoint(std::unique_ptr<geometry::Geometry>&& geometry) override;

        // NOLINTNEXTLINE(modernize-use-nodiscard)
        size_type AddEntity(base::RefEl ref_el,
                            const nonstd::span<const size_type>& nodes,
                            std::unique_ptr<geometry::Geometry>&& geometry) override;
        
        [[nodiscard]] std::shared_ptr<mesh::Mesh> Build() override;

        /** @brief output function printing assembled lists of entity information */
        void PrintLists(std::ostream& o = std::cout) const;

        ~MeshFactory() override = default;

    private:
        dim_t dim_world_;  // dimension of ambient space
        polytopic2d::Mesh::NodeCoordList nodes_;
        polytopic2d::Mesh::EdgeList edges_;
        polytopic2d::Mesh::CellList elements_;

        // If set to true, the Build() method will check whether all sub-entities
        // belong to at least one entity */
        bool check_completeness_;



};

} //namespace lf::mesh::polytoppic2d

#endif // MESH_FACTORY_POLYTOPIC2D_H