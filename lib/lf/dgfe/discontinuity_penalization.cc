/**
 * @file
 * @brief Implementation of discontinutity penalization functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
*/

#include <algorithm>

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

std::shared_ptr<const lf::dgfe::DGFESpace> DiscontinuityPenalization::dgfeSpace(){
    return dgfe_space_ptr_;
}

scalar_t DiscontinuityPenalization::operator()(const lf::mesh::Entity &edge, scalar_t A_f) const {
    LF_VERIFY_MSG(edge.RefEl() == lf::base::RefEl::kSegment(), "Only works for Segments");

    //prepare pointers to adjacent polygons
    auto polygon_pair = dgfe_space_ptr_->AdjacentPolygons(&edge);
    auto polygon_0 = polygon_pair.first.first;
    auto polygon_1 = polygon_pair.second.first;

    //compute areas of biggest simplex in polygon which has "edge" as a Segment
    auto simplex_areas_0 = simplexAreas(*polygon_0, polygon_pair.first.second);
    scalar_t polygon_area_0 = lf::dgfe::integrate(lf::mesh::polytopic2d::Corners(polygon_0), 0, 0);
    auto max_legendre_degree = dgfe_space_ptr_->MaxLegendreDegree();
    //constant appearing in the term p^(2*(d-1)), p is a constant appearing in the definition of shape-regularity => polynomial degrees (one-dimensionel legendre)
    double p_2d_1 = (max_legendre_degree == 1) ? 1 : 4;
    //volume of the edge
    scalar_t edge_volume = lf::geometry::Volume(*(edge.Geometry()));

    if (polygon_1 == nullptr){ //edge is on boundary
        scalar_t max_simplex_area = *std::max_element(simplex_areas_0.begin(), simplex_areas_0.end());

        scalar_t c_inv = c_inv_const_ * std::min(polygon_area_0 / max_simplex_area, p_2d_1);

        // std::cout << "\nBoundary edge. AF = " << A_f << "   p_k_power2 = " << (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree)
        //                 << "       edge volume = " << edge_volume << "     polygon area = " << polygon_area_0 << "       ALles = "
        //                 << A_f * c_sigma_const_ * c_inv * (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree) * edge_volume / polygon_area_0 << "\n";

        return A_f * c_sigma_const_ * c_inv * (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree) * edge_volume / polygon_area_0;

    } else {
        scalar_t max_simplex_area_0 = *std::max_element(simplex_areas_0.begin(), simplex_areas_0.end());

        auto simplex_areas_1 = simplexAreas(*polygon_1, polygon_pair.second.second);
        scalar_t max_simplex_area_1 = *std::max_element(simplex_areas_1.begin(), simplex_areas_1.end());
        scalar_t polygon_area_1 = lf::dgfe::integrate(lf::mesh::polytopic2d::Corners(polygon_1), 0, 0);

        scalar_t c_inv_0 = c_inv_const_ * std::min(polygon_area_0 / max_simplex_area_0, p_2d_1);
        scalar_t c_inv_1 = c_inv_const_ * std::min(polygon_area_1 / max_simplex_area_1, p_2d_1);

        // std::cout << "\nInterior edge. AF = " << A_f << "   p_k_power2 = " << (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree)
        //                 << "       edge volume = " << edge_volume << "     max thingy = " << std::max(c_inv_0 / polygon_area_0, c_inv_1 / polygon_area_1) << "     Alles= "
        //                 << A_f * c_sigma_const_ * (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree) * edge_volume * std::max(c_inv_0 / polygon_area_0, c_inv_1 / polygon_area_1) << "\n";

        return A_f * c_sigma_const_ * (max_legendre_degree + max_legendre_degree) * (max_legendre_degree + max_legendre_degree) * edge_volume * std::max(c_inv_0 / polygon_area_0, c_inv_1 / polygon_area_1);
    }
}

} //namespace lf::dgfe