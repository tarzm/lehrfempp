/**
 * @file
 * @brief Implementation of functionalities for reference basis functions of a discrete space in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#include <lf/fe/fe.h>

#include "legendre_dgfe.h"


namespace lf::dgfe {

scalar_t legendre_polynomial(size_type i, scalar_t x){
    switch(i){
        case 0:
            return 1.0;
        
        case 1:
            return x;
        
        case 2:
            return 1.5 * x * x - 0.5;
        
        default:
            LF_VERIFY_MSG(false, "Only implemented for polynomials of degree 0, 1 and 2");
    }
}

scalar_t legendre_polynomial_dx(size_type n, scalar_t x){
    switch(n){
        case 0:
            return 0.0;
            break;
        
        case 1:
            return 1.0;
            break;
        
        case 2: 
            return 3.0 * x;
            break;
        
        default:
            LF_VERIFY_MSG(false, "Only implemented for polynomials of degree 0, 1 and 2");
    }
}

scalar_t legendre_polynomial_2D(size_type degree_x, size_type degree_y, const Eigen::Vector2d &coord){
    return legendre_polynomial(degree_x, coord[0]) * legendre_polynomial(degree_y, coord[1]);
}

scalar_t legendre_polynomial_2D_dx(size_type degree_x, size_type degree_y, const Eigen::Vector2d &coord){
    switch(degree_x){
        case 0:
            return 0.0;
            break;
        
        case 1:
            return legendre_polynomial(degree_y, coord[1]);
            break;
        
        case 2:
            return 3.0 * coord[0] * legendre_polynomial(degree_y, coord[1]);
            break;
        
        default:
            LF_VERIFY_MSG(false, "Only implemented for polynomials of degree 0, 1 and 2");
    }
}

scalar_t legendre_polynomial_2D_dy(size_type degree_x, size_type degree_y, const Eigen::Vector2d &coord){
    switch(degree_y){
        case 0:
            return 0.0;
            break;
        
        case 1:
            return legendre_polynomial(degree_x, coord[0]);
            break;
        
        case 2:
            return 3.0 * coord[1] * legendre_polynomial(degree_x, coord[0]);
            break;
        
        default:
            LF_VERIFY_MSG(false, "Only implemented for polynomials of degree 0, 1 and 2");
    }
}

scalar_t legendre_basis(size_type n, size_type max_degree, const Eigen::Vector2d &coord){
    auto degrees = multiIndexToDegree(n, max_degree);
    return legendre_polynomial_2D(degrees.first, degrees.second, coord);
}

scalar_t legendre_basis_dx(size_type n, size_type max_degree, const Eigen::Vector2d &coord){
    auto degrees = multiIndexToDegree(n, max_degree);
    return legendre_polynomial_2D_dx(degrees.first, degrees.second, coord);
}

scalar_t legendre_basis_dy(size_type n, size_type max_degree, const Eigen::Vector2d &coord){
    auto degrees = multiIndexToDegree(n, max_degree);
    return legendre_polynomial_2D_dy(degrees.first, degrees.second, coord);
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
