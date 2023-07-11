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
#include <lf/quad/quad.h>
#include <lf/dgfe/dgfe.h>

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
    [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> mesh() { return mesh_; }

    /**
     * @brief Writes the vaues of the MeshFunction in the following format:
     * 
     * x_coordinate y_coordinate value
     * x_coordinate y_coordinate value
     * 
     * points are taken from QR Rule
     */
    template<class MESH_FUNCTION>
    void writeSimple(const MESH_FUNCTION& meshFunction , const std::string& filename){

        std::ofstream out_file;
        out_file.open(filename);

        //quadrule setup
        unsigned integration_degree = 5;
        const lf::quad::QuadRule qr_t = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), integration_degree);
        // qr points
        const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
        //weights
        Eigen::VectorXd w_ref_t{qr_t.Weights()};


        //loop over cells
        for(auto cell : mesh_->Entities(0)){
            //local - global mapping
            lf::dgfe::BoundingBox box(*cell);
            //get sub-tessellation
            auto sub_tessellation = lf::dgfe::subTessellation(cell);

            //loop over triangles in the sub-tessellation
            for(auto& tria_geo_ptr : sub_tessellation){
                //qr points mapped to triangle
                Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
                // qr points mapped back into reference bounding box to retrieve values
                Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
                auto evaluated = meshFunction(*cell, zeta_box_t);

                //loop over points
                for (int i = 0; i < zeta_global_t.cols() ; i++){
                    out_file << zeta_global_t(0,i) << " " << zeta_global_t(1,i) << " " << evaluated[i] << "\n";
                }
            }
        }

        out_file.close();
    }


    /**
     * @brief Create a new NumpyPolytopicWriter. The constructor does not do anything yet
     * @param mesh The mesh
     *
     * @note VtkPolytopicReader supports only writing to ASCII files.
     */
    NumpyPolytopicWriter(std::shared_ptr<const lf::mesh::Mesh> &mesh);


private:

    /// The underlying mesh 
    std::shared_ptr<const lf::mesh::Mesh> mesh_;

};


} // namespace lf::io




#endif //NUMPY_POLYTOPIC_WRITER_H