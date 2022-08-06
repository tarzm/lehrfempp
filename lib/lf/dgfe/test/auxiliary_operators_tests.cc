/**
 * @file
 * @brief Tests of the auxiliary operators functionalities
 * @author Tarzis Maurer
 * @date July 22
 * @copyright ETH Zurich
 */


#include <gtest/gtest.h>

#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::dgfe::test {

TEST(l2ProjectionSqrtANablaBasis, testMesh){

    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    // 2x2 diffusion tensor A(x)
    auto a_simple_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
        return (Eigen::Matrix<double, 2, 2>() << 1.0, 0.0, 0.0, 1.0).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_a_simple{a_simple_lambda};

    auto a_exp_lambda = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
        return (Eigen::Matrix<double, 2, 2>() << exp(x[0] * x[0]), 0.0, 0.0, exp(x[1] * x[1])).finished();
    };
    lf::dgfe::MeshFunctionGlobalDGFE m_a_exp{a_exp_lambda};

    //setup dgfe space
    unsigned max_legendre_degree = 2;
    std::shared_ptr<lf::dgfe::DGFESpace> dgfe_space_ptr(new lf::dgfe::DGFESpace(mesh_ptr, max_legendre_degree));

    auto l2_projection_simple = L2ProjectionSqrtANablaBasis<double>(dgfe_space_ptr, m_a_simple, 20);
    auto l2_projection_exp = L2ProjectionSqrtANablaBasis<double>(dgfe_space_ptr, m_a_exp, 20);

    //!!!!!!!! SETUP MESH FUNCTIONS!!!!!!!!!!!!
    unsigned basis;
    //mesh function lambda for calculated projection
    auto m_l2_projection_simple = [&m_a_simple, &basis, &l2_projection_simple](const lf::mesh::Entity &cell, Eigen::MatrixXd local) -> std::vector<Eigen::Vector2d> {
        std::vector<Eigen::Vector2d> result;

        for (int i = 0; i < local.cols(); i++){
            Eigen::Vector2d res = { l2_projection_simple.first[basis](cell, local.col(i))[0],
                                    l2_projection_simple.second[basis](cell, local.col(i))[0]};
            result.push_back(res);
        }
        return result;
    };
    //mesh function lambda for true function
    auto true_simple = [&m_a_simple, &basis, &max_legendre_degree](const lf::mesh::Entity &cell, Eigen::MatrixXd local) -> std::vector<Eigen::Vector2d> {
        std::vector<Eigen::Vector2d> result;
        auto a_eval = m_a_simple(cell, local);
        lf::dgfe::BoundingBox box(cell);
        
        for (int i = 0; i < local.cols(); i++){
            Eigen::Vector2d nabla_basis{legendre_basis_dx(basis, max_legendre_degree, local.col(i)) * box.inverseJacobi(0),
                                        legendre_basis_dy(basis, max_legendre_degree, local.col(i)) * box.inverseJacobi(1)};
            result.push_back(a_eval[i].sqrt() * nabla_basis);
        }

        return result;
    };

    //mesh function lambda for calculated projection
    auto m_l2_projection_exp = [&m_a_exp, &basis, &l2_projection_exp](const lf::mesh::Entity &cell, Eigen::MatrixXd local) -> std::vector<Eigen::Vector2d> {
        std::vector<Eigen::Vector2d> result;

        for (int i = 0; i < local.cols(); i++){
            Eigen::Vector2d res = { l2_projection_exp.first[basis](cell, local.col(i))[0],
                                    l2_projection_exp.second[basis](cell, local.col(i))[0]};
            result.push_back(res);
        }
        return result;
    };
    //mesh function lambda for true function
    auto true_exp = [&m_a_exp, &basis, &max_legendre_degree](const lf::mesh::Entity &cell, Eigen::MatrixXd local) -> std::vector<Eigen::Vector2d> {
        std::vector<Eigen::Vector2d> result;
        auto a_eval = m_a_exp(cell, local);
        lf::dgfe::BoundingBox box(cell);
        
        for (int i = 0; i < local.cols(); i++){
            Eigen::Vector2d nabla_basis{legendre_basis_dx(basis, max_legendre_degree, local.col(i)) * box.inverseJacobi(0),
                                        legendre_basis_dy(basis, max_legendre_degree, local.col(i)) * box.inverseJacobi(1)};
            result.push_back(a_eval[i].sqrt() * nabla_basis);
        }

        return result;
    };
    //!!!!!!!! END SETUP MESH FUNCTIONS!!!!!!!!!!!!

    //calculate errors for each basis function
    unsigned n_local_dofs = (max_legendre_degree + 1) * (max_legendre_degree + 1);

    for (basis = 1; basis < n_local_dofs; basis++){
        double l2_error_simple = L2ErrorGradSubTessellation<double>(m_l2_projection_simple, true_simple, mesh_ptr, 40);
        double l2_error_exp = L2ErrorGradSubTessellation<double>(m_l2_projection_exp, true_exp, mesh_ptr, 40);

        EXPECT_NEAR(l2_error_simple, 0.0, std::numeric_limits<double>::epsilon() * 100) << " At " << basis << "\n";
        EXPECT_TRUE(l2_error_simple * 100000 < l2_error_exp) << " At " << basis << "\n";
        EXPECT_TRUE(l2_error_exp < 0.1) << " At " << basis << "\n";

        std::cout << "L2 Error of exponential diffusion coefficient is " << l2_error_exp << "\n";
    }
}

}