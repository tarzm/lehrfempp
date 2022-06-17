/**
 * @file
 * @brief Implementation sub-tessellation functionalities
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#include "sub_tessellation.h"

namespace lf::dgfe {

using scalar_t = double;

std::vector<std::unique_ptr<lf::geometry::TriaO1>> subTessellation(const lf::mesh::Entity *polygon){
    LF_VERIFY_MSG(polygon->RefEl() == lf::base::RefEl::kPolygon(), "Sub-tessellation only works with Polygons");

    std::vector<std::unique_ptr<lf::geometry::TriaO1>> result;
    
    //get cells nodes
    auto corners = lf::mesh::polytopic2d::Corners(polygon);
    //get barycenter of cell
    Eigen::MatrixXd barycenter = corners.rowwise().mean();

    //number of nodes and segments is equal in polygons
    result.reserve(corners.cols());

    //loop over segments of the polygon
    for (auto segment : polygon->SubEntities(1)){
        //coordinate matrix for triangle to be generated
        Eigen::MatrixXd tria_corners(2,3);
        //coordinates of the segments endpoints
        Eigen::MatrixXd segment_corners = lf::geometry::Corners(*(segment->Geometry()));
        tria_corners.block(0,0,2,2) = segment_corners;
        tria_corners.col(2) = barycenter.col(0);
        result.emplace_back(std::make_unique<lf::geometry::TriaO1>(tria_corners));
    }

    return result;
}


} //namespace lf::dgfe