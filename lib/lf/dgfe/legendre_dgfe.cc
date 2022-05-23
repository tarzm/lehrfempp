/**
 * @file
 * @brief Implementation of classes supporting the discontinuous finite element discretization of boundary value problems
 * with mapped legendre basis functions
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/



#include "legendre_dgfe.h"


namespace lf::dgfe {

scalar_t legendre_polynomial(size_type i, scalar_t x){
    return legendre_coeffs_(i,0) + legendre_coeffs_(i,1) * x + legendre_coeffs_(i,2) * std::pow(x, 2);
}

scalar_t C_i_j_k(size_type i, size_type j, size_type k, size_type degree_p){
    scalar_t sum = 0;
    for (int n = 0; n <= degree_p; n++){
        for (int m = 0; m <= degree_p; m++){
            if( n+m == k){
                sum += legendre_coeffs_(i, n) * legendre_coeffs_(j, m);
            }
        }
    }
    return sum;
}

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





} //namespace lf::dgfe
