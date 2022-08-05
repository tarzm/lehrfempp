/**
 * @file
 * @brief Declares the NumpyPolytopicWriter which can read vtk files defining polytopic meshes
 * @author Tarzis Maurer
 * @date   2022-05-13
 * @copyright ETH Zurich
 */

#ifndef NUMPY_POLYTOPIC_WRITER_H
#define NUMPY_POLYTOPIC_WRITER_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>
#include <fstream>

namespace lf::io {

using scalar_t = double;

class NumpyPolytopicWriter {

public:

    using size_type = lf::mesh::Mesh::size_type;
    using dim_t = lf::base::RefEl::dim_t;

    /**
     * @brief Get the mesh that is to be written
     */
    [[nodiscard]] std::shared_ptr<mesh::Mesh> mesh() { return mesh_; }

    /**
     * @brief Writes the mesh in the following format:
     * 
     * x_coordinate y_coordinate value
     * x_coordinate y_coordinate value
     * 
     * At first for all nodes, then for all baricenters of all cells
     */
    template<class MESH_FUNCTION>
    void writeSimple(const MESH_FUNCTION& meshFunction , const std::string& filename){

        std::ofstream out_file;
        out_file.open(filename);

        //all nodes
        for (auto point : mesh_->Entities(2)){
            auto corners = lf::geometry::Corners(*(point->Geometry()));
            Eigen::Matrix<double, 0, 1> dummy;
            auto value_vec = meshFunction(*point, dummy); //this is a vector
            out_file << corners(0,0) << " " << corners(1,0) << " " << value_vec[0] << "\n";
        }  
        //baricenters of all cells
        for(auto cell : mesh_->Entities(0)){
            auto corners = lf::mesh::polytopic2d::Corners(cell);
            Eigen::MatrixXd mean = corners.rowwise().mean();
            auto value_vec = meshFunction(*cell, mean); //this is a vector
            out_file << mean(0,0) << " " << mean(1,0) << " " << value_vec[0] << "\n";
        }




        out_file.close();
    }


    /**
     * @brief Create a new NumpyPolytopicWriter. The constructor does not do anything yet
     * @param mesh The mesh
     *
     * @note VtkPolytopicReader supports only writing to ASCII files.
     */
    NumpyPolytopicWriter(std::shared_ptr<lf::mesh::Mesh> &mesh);


private:

    /// The underlying mesh 
    std::shared_ptr<mesh::Mesh> mesh_;

};


} // namespace lf::io




#endif //NUMPY_POLYTOPIC_WRITER_H