/**
 * @file
 * @brief Implementation of discontinutity penalization functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
*/

#include "discontinuity_penalization.h"

namespace lf::dgfe {

std::vector<scalar_t> simplexAreas(const lf::mesh::Entity &cell, size_type egde_loc_idx){
    LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kPolygon(), "Only implemented for polygons");
    
    std::vector<scalar_t> result;
    auto corners = lf::mesh::polytopic2d::Corners(&cell);
    auto n_nodes = corners.cols();

    //get local indices of two endpoints
    size_type idx_0 = egde_loc_idx;
    size_type idx_1 = (egde_loc_idx + 1) % n_nodes;

    Eigen::MatrixXd coords = Eigen::MatrixXd::Zero(2,3);
    coords.col(0) = corners.col(idx_0);
    coords.col(1) = corners.col(idx_1);

    for (int i = 1; i < n_nodes - 1; i++){
        size_type next_node_idx = (idx_1 + i) % n_nodes;
        coords.col(2) = corners.col(next_node_idx);
        result.push_back(lf::dgfe::integrate(coords, 0, 0));
    }

    LF_VERIFY_MSG(result.size() == n_nodes - 2, "Something went wrong with the running variable i");
    return result;
}








} //namespace lf::dgfe