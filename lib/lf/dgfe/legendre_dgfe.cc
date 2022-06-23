/**
 * @file
 * @brief Implementation of functionalities for reference basis functions of a discrete space in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/



#include "legendre_dgfe.h"


namespace lf::dgfe {

scalar_t legendre_polynomial(size_type i, scalar_t x){
    return legendre_coeffs_(i,0) + legendre_coeffs_(i,1) * x + legendre_coeffs_(i,2) * x * x;
}

scalar_t legendre_polynomial_2D(size_type degree_x, size_type degree_y, Eigen::Vector2d coord){
    return legendre_polynomial(degree_x, coord[0]) * legendre_polynomial(degree_y, coord[1]);
}


scalar_t C_i_j_k(size_type i, size_type j, size_type k){
    scalar_t sum = 0;
    for (int n = 0; n <= i; n++){
        for (int m = 0; m <= j; m++){
            if( n+m == k){
                sum += legendre_coeffs_(i, n) * legendre_coeffs_(j, m);
            }
        }
    }
    return sum;
}

std::pair<size_type, size_type> multiIndexToDegree(size_type basis_index, size_type degree_of_space){
    size_type degree_x;
    size_type degree_y;
    switch(degree_of_space){
        case 1:
            LF_VERIFY_MSG(basis_index < 4 && basis_index >= 0, "Basis index must be 0 <= basis_index <= 4 when degree of the space is 1");
            degree_x = basis_index / 2;
            degree_y = basis_index % 2;
            return std::make_pair(degree_x, degree_y);

        case 2:
            LF_VERIFY_MSG(basis_index < 9 && basis_index >= 0, "Basis index must be 0 <= basis_index <= 8 when degree of the space is 2");
            degree_x = basis_index / 3;
            degree_y = basis_index % 3;
            return std::make_pair(degree_x, degree_y);
        default:
            LF_VERIFY_MSG(false, "Degree of space must be 0 <= degree <= 2");
    }
}

size_type degreeToMultiIndex(std::pair<size_type, size_type> degrees, size_type degree_of_space){
    switch(degree_of_space){
        case 1:
            LF_VERIFY_MSG(degrees.first < 2 && degrees.first >= 0 && degrees.second < 2 && degrees.second >= 0, "Degrees must be 0 <= degree <= 1 when space is of degree 1");
            return degrees.second + degrees.first * 2;
        case 2:
            LF_VERIFY_MSG(degrees.first < 3 && degrees.first >= 0 && degrees.second < 3 && degrees.second >= 0, "Degrees must be 0 <= degree <= 2 when space is of degree 2");
            return degrees.second + degrees.first * 3;
        default:
            LF_VERIFY_MSG(false, "Degree of space must be 0 <= degree <= 2");
    }
}




} //namespace lf::dgfe
