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
            auto degrees_i = multiIndexToDegree(i, 1);
            auto degrees_j = multiIndexToDegree(j, 1);
            int i1 = degrees_i.first;   //degree of x in trial basis
            int i2 = degrees_i.second;  //degree of y in trial basis
            int j1 = degrees_j.first;   //degree of x in test basis
            int j2 = degrees_j.second;  //degree of y in test basis

            for (int k = 0; k <= i1 + j1; k++){
                for (int l = 0; l <= i2 + j2; l++){
                    scalar_t C_i1_j1_k = C_i_j_k(i1, j1, k);
                    scalar_t C_i2_j2_l = C_i_j_k(i2, j2, l);
                    if (C_i1_j1_k > 0 && C_i2_j2_l > 0){
                        sum += C_i1_j1_k * C_i2_j2_l * box.det() * lf::dgfe::integrate(box.map(corners), k, l);
                    }
                    
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
            auto degrees = multiIndexToDegree(basis, 1);
            int legendre_degree_x = degrees.first;
            int legendre_degree_y = degrees.second;
            //legendre polynomials have max degree 1, therefore no legendre_coefficients have to be taken into account
            elem_vec[basis] += monomial_coeff * box.det() * lf::dgfe::integrate(box.map(corners), legendre_degree_x + monomial_degree_x, legendre_degree_y + monomial_degree_y);

            //std::cout << "Added term " << monomial_coeff * box.det() * lf::dgfe::integrate(box.map(corners), legendre_degree_x + monomial_degree_x, legendre_degree_y + monomial_degree_y) << "\n";
        }
    }

    return elem_vec;
}

Eigen::Matrix<scalar_t, 9, 9> DGFEO2MassElementMatrix::Eval(const lf::mesh::Entity &cell) const{
    BoundingBox box(cell);
    auto corners = lf::mesh::polytopic2d::Corners(&cell);
    Eigen::Matrix<scalar_t, 9, 9> elem_mat;
    
    //loop over trial basis funtions on cell
    for (int i = 0; i < 9; i++){
        //loop over test basis functions on cell
        for (int j = 0; j < 9; j++){
            scalar_t sum = 0;
            //definition of i1, i2, j1, j2
            auto degrees_i = multiIndexToDegree(i, 2);
            auto degrees_j = multiIndexToDegree(j, 2);
            int i1 = degrees_i.first;   //degree of x in trial basis
            int i2 = degrees_i.second;  //degree of y in trial basis
            int j1 = degrees_j.first;   //degree of x in test basis
            int j2 = degrees_j.second;  //degree of y in test basis

            for (int k = 0; k <= i1 + j1; k++){
                for (int l = 0; l <= i2 + j2; l++){
                    scalar_t C_i1_j1_k = C_i_j_k(i1, j1, k);
                    scalar_t C_i2_j2_l = C_i_j_k(i2, j2, l);
                    if (C_i1_j1_k > 0 && C_i2_j2_l > 0){
                        sum += C_i1_j1_k * C_i2_j2_l * box.det() * lf::dgfe::integrate(box.map(corners), k, l);
                    }
                    
                }
            }
            elem_mat(i,j) = sum;
        }
    }
    return elem_mat;
}

DGFEO2LocalLoadVector::ElemVec DGFEO2LocalLoadVector::Eval(const lf::mesh::Entity &cell) const {
    BoundingBox box(cell);
    auto corners = lf::mesh::polytopic2d::Corners(&cell);

    ElemVec elem_vec = Eigen::Matrix<scalar_t, 9, 1>::Zero();

    //loop over monomials of the polynomial
    for (auto monomial : polynomial_){
        auto monomial_coeff = monomial.first;
        auto monomial_degree_x = monomial.second.first;
        auto monomial_degree_y = monomial.second.second;

        //loop over trial basis funtions on cell
        for(int basis = 0; basis < 9; basis++){
            auto basis_degrees = multiIndexToDegree(basis, 2);
            int legendre_degree_x = basis_degrees.first;
            int legendre_degree_y = basis_degrees.second;
        
            for (int degree_x = 0; degree_x <= legendre_degree_x; degree_x++){
                for (int degree_y = 0; degree_y <= legendre_degree_y; degree_y++){
                    //check wether coefficients are not 0
                    scalar_t coeff_x = legendre_coeffs_(legendre_degree_x, degree_x);
                    scalar_t coeff_y = legendre_coeffs_(legendre_degree_y, degree_y);
                    if (coeff_x != 0.0 && coeff_y != 0.0){
                        elem_vec[basis] += monomial_coeff * box.det() * coeff_x * coeff_y *
                                            lf::dgfe::integrate(box.map(corners), degree_x + monomial_degree_x, degree_y + monomial_degree_y);
                    }
                }
            }
        }
    }

    return elem_vec;
}

DGFEMassElementMatrixST::DGFEMassElementMatrixST(unsigned max_integration_degree, unsigned max_legendre_degree) : max_integration_degree_(max_integration_degree), max_legendre_degree_(max_legendre_degree) {
    LF_VERIFY_MSG(max_legendre_degree_ == 1 || max_legendre_degree_ == 2, "Only implemented for maximum 1D legendre polynomials of degree 1 and 2");
}

Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> DGFEMassElementMatrixST::Eval(const lf::mesh::Entity &cell) const {

    unsigned matrix_size = (max_legendre_degree_ == 1) ? 4 : 9;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> elem_mat(matrix_size, matrix_size);

    auto corners = lf::mesh::polytopic2d::Corners(&cell);

    int i1;
    int i2;
    int j1;
    int j2;
    auto eval_lambda = [&i1, &i2, &j1, &j2](const lf::mesh::Entity *entity, Eigen::Vector2d coord) -> scalar_t{
        return lf::dgfe::legendre_polynomial_2D(i1, i2, coord) * lf::dgfe::legendre_polynomial_2D(j1, j2, coord);
    };

    lf::dgfe::SubTessellationIntegrator<scalar_t, decltype(eval_lambda)> integrator;
    

    //loop over trial basis funtions on cell
    for (int i = 0; i < matrix_size; i++){
        //loop over test basis functions on cell
        for (int j = 0; j < matrix_size; j++){
            scalar_t sum = 0;
            //definition of i1, i2, j1, j2
            auto degrees_i = multiIndexToDegree(i, max_legendre_degree_);
            auto degrees_j = multiIndexToDegree(j, max_legendre_degree_);
            i1 = degrees_i.first;   //degree of x in trial basis
            i2 = degrees_i.second;  //degree of y in trial basis
            j1 = degrees_j.first;   //degree of x in test basis
            j2 = degrees_j.second;  //degree of y in test basis

            elem_mat(i,j) = integrator.integrate(&cell, eval_lambda, max_integration_degree_);
        }
    }

    return elem_mat;
}

} //namespace lf::dgfe