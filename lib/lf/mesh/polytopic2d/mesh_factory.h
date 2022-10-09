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
         * @param unit_square If set to true, the MeshFactory will "clean" the mesh from coordinates that are very close to 1 and 0
         * but not exactly. It will make them exactly 0 and 1. This problem is coming from PolyMesher, which collepses small edges
         * into a single point. After that, the polygons do not make up exactly a unit square anymore but form a slightly curved boundary,
         * which causes errors of integration over the mesh to rise significantly (up to 7 orders of magnitude observed for a polynomial integrand).
         */
        explicit MeshFactory(dim_t dim_world, bool check_completeness = true, bool unit_square = true)
            : dim_world_(dim_world), check_completeness_(check_completeness), unit_square_(unit_square) {}

        
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
        polytopic2d::Mesh::CellList elements_;

        // If set to true, the Build() method will check whether all sub-entities
        // belong to at least one entity */
        bool check_completeness_;
        bool unit_square_;



};

/**
 * @brief returns a polytopic 2D mesh from a hybrid 2D mesh
 * 
 * @return lf::mesh::polytopic2d::Mesh 
 */
std::shared_ptr<lf::mesh::Mesh> polytopicFromHybrid2D(std::shared_ptr<const lf::mesh::Mesh> mesh_ptr);

} //namespace lf::mesh::polytoppic2d

#endif // MESH_FACTORY_POLYTOPIC2D_H