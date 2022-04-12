/**
 * @file
 * @brief Declares the VtkPolytopicReader which can read vtk files defining polytopic meshes
 * @author Tarzis Maurer
 * @date   2022-04-12
 * @copyright ETH Zurich
 */

#ifndef VTK_POLYTOPIC_READER_H
#define VTK_POLYTOPIC_READER_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/utils/utils.h>


#include <string>

namespace lf::io {

enum VtkSection{
    header, //0
    points, //1
    metadata, //2
    polygons, //3
    offsets, //4
    connectivity, //5
};

class VtkPolytopicReader {

    public:

        using size_type = lf::mesh::Mesh::size_type;
        using dim_t = lf::base::RefEl::dim_t;

        /**
         * @brief Get the mesh that was read by this reader.
         */
        [[nodiscard]] std::shared_ptr<mesh::Mesh> mesh() { return mesh_; }

        /**
         * @brief Get the mesh that was read by this reader.
         */
        [[nodiscard]] std::shared_ptr<const mesh::Mesh> mesh() const { return mesh_; }

        /**
         * @brief Create a new VtkPolytopicReader by reading from the specified file.
         * @param factory The mesh::MeshFactory that is used to construct the mesh.
         * @param filename The filename of the `.vtk` file that is read.
         *
         * @note VtkPolytopicReader supports only ASCII `.vtk` files.
         */
        VtkPolytopicReader(std::unique_ptr<mesh::polytopic2d::MeshFactory> factory,
                    const std::string& filename);



    private:

        /// The underlying mesh created by the mesh factory.
        std::shared_ptr<mesh::Mesh> mesh_;

        //mesh factory to create the mesh
        std::unique_ptr<mesh::MeshFactory> mesh_factory_;

};









} //namespace lf::io

#endif // VTK_POLYTOPICV_READER_H