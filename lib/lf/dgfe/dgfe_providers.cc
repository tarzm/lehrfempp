/**
 * @file
 * @brief Implementation of classes supporting assembly of element matrices/vectors in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#include "dgfe_providers.h"

namespace lf::dgfe {

Eigen::Matrix<scalar_t, 4, 4> DGFEO1MassElementMatrix::Eval(const lf::mesh::Entity &cell) const{
    BoundingBox box(cell);
    auto corners = lf::mesh::polytopic2d::Corners(&cell);
    Eigen::Matrix<scalar_t, 4, 4> elem_mat;
    
    //loop over trial basis funtions on cell
    for (int i = 0; i < 4; i++){
        //loop over test basis functions on cell
        for (int j = 0; j < 4; j++){
            scalar_t sum = 0;
            //definition of i1, i2, j1, j2
            int i1 = i / 2;
            int i2 = i % 2;
            int j1 = j / 2;
            int j2 = j % 2;
            for (int k = 0; k <= i1 + j1; k++){
                for (int l = 0; l <= i2 + j2; l++){
                    sum += C_i_j_k(i1, j1, k, 1) * C_i_j_k(i2, j2, l, 1) * box.det() * lf::dgfe::integrate(box.map(corners), k, l);
                }
            }
            elem_mat(i,j) = sum;
        }
    }
    return elem_mat;
}

DGFEO1LocalLoadVector::ElemVec DGFEO1LocalLoadVector::Eval(const lf::mesh::Entity &cell) const {
    BoundingBox box(cell);
    auto corners = lf::mesh::polytopic2d::Corners(&cell);

    ElemVec elem_vec = Eigen::Matrix<scalar_t, 4, 1>::Zero();

    //loop over monomials of the polynomial
    for (auto monomial : polynomial_){
        auto monomial_coeff = monomial.first;
        auto monomial_degree_x = monomial.second.first;
        auto monomial_degree_y = monomial.second.second;

        //loop over trial basis funtions on cell
        for(int basis = 0; basis < 4; basis++){
            int basis_degree_x = basis / 2;
            int basis_degree_y = basis % 2;
            //legendre polynomials have max degree 1, therefore no legendre_coefficients have to be taken into account
            elem_vec[basis] += monomial_coeff * box.det() * lf::dgfe::integrate(box.map(corners), basis_degree_x + monomial_degree_x, basis_degree_y + monomial_degree_y);

            std::cout << "Added term " << monomial_coeff * box.det() * lf::dgfe::integrate(box.map(corners), basis_degree_x + monomial_degree_x, basis_degree_y + monomial_degree_y) << "\n";
        }
    }

    return elem_vec;
}

} //namespace lf::dgfe