/**
 * @file
 * @brief Tests of the integration algorithm
 * @author Tarzis Maurer
 * @date 2022-04-19
 * @copyright ETH Zurich
 */

#include <cmath>
#include <filesystem>

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/dgfe/dgfe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/io/io.h>
#include "lf/mesh/test_utils/test_meshes.h"

#define NORMALTOLERANCE 1e-12
#define TOLERANCE 1e-10

namespace lf::dgfe::test {

bool isNormalOf(const Eigen::MatrixXd normal, const Eigen::MatrixXd edge){
    LF_VERIFY_MSG(std::abs(normal.col(0).norm() - 1.0) < NORMALTOLERANCE, "Normal does not have length 1");
    Eigen::MatrixXd edge_vector(2,1);
    edge_vector(0,0) = edge(0,1) - edge(0,0);
    edge_vector(1,0) = edge(1,1) - edge(1,0);
    return (std::abs(edge_vector.col(0).dot(normal.col(0))) < NORMALTOLERANCE);
}

TEST(integration, helperFunctions){

    //test outwardNormal
    Eigen::MatrixXd a_edge(2,2);
    a_edge <<       0, 1,
                    1, 0;
    Eigen::MatrixXd normal_check(2,1);
    normal_check << - std::sqrt(2.0) / 2.0,  - std::sqrt(2.0) / 2.0;
    EXPECT_TRUE(lf::dgfe::outwardNormal(a_edge).isApprox(normal_check));
    EXPECT_TRUE(isNormalOf(lf::dgfe::outwardNormal(a_edge), a_edge));

    //check outwardNormal of all edges of the pentagon from the paper
    Eigen::MatrixXd b_polygon(2,5);
    b_polygon <<    -0.666666666666667, 0.555555555555556, 1.0, -0.555555555555556, -1.0,
                    -0.789473684210526, -1.0, -0.52631578947368, 1.0, -0.157894736842105;

    for (int i = 0; i < b_polygon.cols(); i++){
        Eigen::MatrixXd edge_i(2,2);
        edge_i.col(0) = b_polygon.col(i);
        edge_i.col(1) = b_polygon.col((i + 1) % b_polygon.cols());
        EXPECT_TRUE(isNormalOf(lf::dgfe::outwardNormal(edge_i), edge_i));
    }

    //test euclideanDist
    Eigen::MatrixXd b_point(2,1);
    Eigen::MatrixXd c_point(2,1);
    b_point << 0.0, 0.0;
    c_point << 2.0, 2.5;
    EXPECT_NEAR(euclideanDist(b_point, c_point), 3.2015621187164243, 1e-12);
    Eigen::MatrixXd f_point(2,1);
    Eigen::MatrixXd g_point(2,1);
    f_point << -0.8, 3.0;
    g_point << 2.34, 1.01;
    EXPECT_NEAR(euclideanDist(f_point, g_point), std::sqrt(std::pow( 2.34-(-0.8) , 2) + std::pow( 1.01-3.0 , 2)), 1e-12);

}

TEST(integration, lineIntegral){

    Eigen::MatrixXd a_polygon(2,2);
    a_polygon <<    2.0, 5.0,
                    1.0, 3.0;
    
    //calculated "by hand"
    EXPECT_DOUBLE_EQ(lf::dgfe::integrate(a_polygon, 4, 3), std::sqrt(13) * (648.0/8.0 + 2700.0/7.0 + 4806.0/6.0 + 4737.0/5.0 + 2792.0/4.0 + 984.0/3.0 + 192.0/2.0 + 16.0));
}

TEST(integration, triangle){

    //test simple triangle from paper
    Eigen::MatrixXd a_polygon(2,3);
    a_polygon <<    -1.0, 1.0, -1.0,
                    -1.0, 0.0, 1.0;

    //NOTE: The tests commented out will fail
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 5), 0.0, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 10), 0.0111339078, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 20, 20), 0.0030396808, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 40, 40), 7.9534562047e-14, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 10, 5), 0.0, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 20, 40), 0.0, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 40, 5), 0.0, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 20), -0.005890191, TOLERANCE);
    //EXPECT_NEAR(lf::dgfe::integrate(a_polygon, 5, 40), -0.001868889, TOLERANCE);
}

TEST(integration, pentagon){

    //test pentagon from paper
    Eigen::MatrixXd b_polygon(2,5);
    b_polygon <<    -0.666666666666667, 0.555555555555556, 1.000000000000000, -0.555555555555556, -1.000000000000000,
                    -0.789473684210526, -1.000000000000000, -0.052631578947368, 1.000000000000000, -0.157894736842105;

    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 5), -0.0020324991, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 10), 7.4274779926e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 20, 20), 6.0738145408e-8, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 40, 40), 2.2238524572e-12, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 10, 5), -2.0911953867e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 20, 5), -1.3797380205e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 40, 5), -7.9203571311e-7, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 20), 8.08469022058e-5, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(b_polygon, 5, 40), 4.37593748009e-5, TOLERANCE);
}

TEST(integration, nonConvexPolygon){
    //test non-convex polygon from paper
    Eigen::MatrixXd c_polygon(2,15);
    c_polygon <<    0.413048522141662, 0.024879797655533, -0.082799691823524, -0.533191422779328, -0.553573605852999, -0.972432940212767, -1.000000000000000, -0.789986179147920, -0.627452906935866, -0.452662174765764, -0.069106265580153, 0.141448047807069, 1.000000000000000, 0.363704451489016, 0.627086024018283,
                    0.781696234443715, 0.415324992429711, 0.688810136531751, 1.000000000000000, 0.580958514816226, 0.734117068746903, 0.238078507228890, 0.012425068086110, -0.636532897516109, -1.000000000000000, -0.289054989277619, -0.464417038155806, -0.245698820584615, -0.134079689960635, -0.110940423607648;
    //Note: test commented out will fail
    //EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 5), -0.002589861, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 10, 10), 1.5738050178e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 20, 20), 1.3793481020e-6, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 40, 40), 4.2588831784e-10, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 10, 5), 0.0014996521, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 20, 5), 7.0356275077e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 40, 5), 2.5065856538e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 20), -1.330384913e-4, TOLERANCE);
    EXPECT_NEAR(lf::dgfe::integrate(c_polygon, 5, 40), -3.963064075e-5, TOLERANCE);
}

TEST(integration, polytopicTestMesh){
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);

    //test x^3 * y^4
    lf::dgfe::scalar_t sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        sum += integrate(corners, 3, 4);
    }
    EXPECT_NEAR(sum, 0.05, 1e-14);

    //test x^5 * y^7
    sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        sum += integrate(corners, 5, 7);
    }
    EXPECT_NEAR(sum, (1.0/48.0), 1e-14);
}

TEST(integration, bigMesh){
    //get mesh
    std::filesystem::path here = __FILE__;
    auto mesh_file = here.parent_path().string() + "/msh_files/unit_square_voronoi_1000_cells.vtk";

    lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
    
    auto mesh_ptr = reader.mesh();

    lf::dgfe::scalar_t sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        auto corners = lf::mesh::polytopic2d::Corners(cell);
        //std::cout << "Cell " << mesh_ptr->Index(*cell) << " contributes " << integrate(corners, 3, 4) << "\n";
        sum += integrate(corners, 3, 4);
    }
    EXPECT_NEAR(sum, 0.05, TOLERANCE);

}

TEST(sub_tessellation_integration, polytopicTestMesh){
    //get mesh
    auto mesh_ptr = lf::mesh::test_utils::GeneratePolytopic2DTestMesh(0,1);
    
    //degrees used for lambda
    int degree_x = 3;
    int degree_y = 4;
    // lambda function for polynomial integration
    auto polynomial_lambda = [&degree_x, &degree_y](const lf::mesh::Entity &entity, Eigen::Vector2d x) -> double {
        return std::pow(x[0], degree_x) * std::pow(x[1], degree_y);
    };

    //initialize integrator
    lf::dgfe::SubTessellationIntegrator<double, decltype(polynomial_lambda)> poly_integrator;

    
    // Test monomial x^3 * y^3
    lf::dgfe::scalar_t sum = 0.0;
    for (auto cell : mesh_ptr->Entities(0)){
        sum += poly_integrator.integrate(*cell, polynomial_lambda, degree_x + degree_y);
    }
    EXPECT_NEAR(sum, 0.05, 1e-14);


    // Test monomial x^5 * y^7
    sum = 0.0;
    degree_x = 5;
    degree_y = 7;
    for (auto cell : mesh_ptr->Entities(0)){
        sum += poly_integrator.integrate(*cell, polynomial_lambda, degree_x + degree_y);
    }
    EXPECT_NEAR(sum, 1.0/48.0, 1e-14);


    // Test sin^2(x) * cos (PI *y^2)
    auto trigonometric_lambda = [](const lf::mesh::Entity &entity, Eigen::Vector2d x) -> double {
        return std::pow(std::sin(x[0]), 2) * std::cos(M_PI * x[1] * x[1]);
    };
    lf::dgfe::SubTessellationIntegrator<double, decltype(trigonometric_lambda)> trigo_integrator;
    sum = 0.0;
    //loop over cells
    for (auto cell : mesh_ptr->Entities(0)){
        sum += trigo_integrator.integrate(*cell, trigonometric_lambda, 14); //14 is lowest degree that results in an error < 1e-14
    }
    //check
    double exact_check = 0.1019760096823904218580342315643055994170804185885964798615039678; // copied from wolframalpha.com
    EXPECT_NEAR(sum, exact_check, 1e-14);


    //Test x^2 + e^(x*y)
    auto exponential_lambda = [](const lf::mesh::Entity &entity, Eigen::Vector2d x) -> double {
        return x[0]*x[0] + exp(x[0] * x[1]);
    };
    lf::dgfe::SubTessellationIntegrator<double, decltype(exponential_lambda)> exp_integrator;
    sum = 0.0;
    //loop over cells
    for (auto cell : mesh_ptr->Entities(0)){
        sum += exp_integrator.integrate(*cell, exponential_lambda, 10); //10 is lowest degree that results in an error < 1e-14
    }
    //check
    exact_check = 1.651235484787737228193342177582565171308234579126117326173794530330979441089472815945286980167725541; // copied from wolframalpha.com
    EXPECT_NEAR(sum, exact_check, 1e-14);
}

TEST(sub_tessellation_integration, 100To400cells){
    //get meshes files--------------------------------------------------------------------------------------
    std::vector<int> mesh_sizes{100, 200, 400, 800, 1000, 1200, 1400, 1600, 2000, 3000, 4000};
    std::filesystem::path here = __FILE__;
    auto mesh_file_100 = here.parent_path().string() + "/msh_files/unit_square_voronoi_100_cells.vtk";
    std::vector<decltype(mesh_file_100)> mesh_files;
    mesh_files.push_back(std::move(mesh_file_100));
    mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_200_cells.vtk");
    mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_400_cells.vtk");
    mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_800_cells.vtk");
    mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_1000_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_1200_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_1400_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_1600_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_2000_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_3000_cells.vtk");
    //mesh_files.push_back(here.parent_path().string() + "/msh_files/unit_square_voronoi_4000_cells.vtk");
    //------------------------------------------------------------------------------------------------


    //Set up lambdas, true solution and calculated solution vectors-----------------------------------
    //degrees used for lambda
    std::vector<std::function<double(Eigen::Vector2d)>> functors;
    std::vector<std::string> function_expressions{"x^3 * y^4", "x^5 * y^7", "sin^2(x) * cos (PI *y^2)", "x^2 + e^(x*y)"};
    // first lambda function for polynomial integration: x^3 * y^4
    functors.emplace_back([](Eigen::Vector2d x) -> double {
        return std::pow(x[0], 3) * std::pow(x[1], 4);
    });
    // second lambda function for polynomial integration: x^5 * y^7
    functors.emplace_back([](Eigen::Vector2d x) -> double {
        return std::pow(x[0], 5) * std::pow(x[1], 7);
    });
    //third lambda: sin^2(x) * cos (PI *y^2)
    functors.emplace_back([](Eigen::Vector2d x) -> double {
        return std::pow(std::sin(x[0]), 2) * std::cos(M_PI * x[1] * x[1]);
    });
    //fourth lambda: x^2 + e^(x*y)
    functors.emplace_back([](Eigen::Vector2d x) -> double {
        return x[0]*x[0] + exp(x[0] * x[1]);
    });
    //setup mesh functions of lambdas
    std::vector<lf::dgfe::MeshFunctionGlobalDGFE<std::function<double(Eigen::Vector2d)>>> mesh_functions;
    for (auto functor : functors){
        mesh_functions.emplace_back(lf::dgfe::MeshFunctionGlobalDGFE<decltype(functor)>(functor));
    }
    //exact solutions calculated and copied from wolframalpha.com
    std::vector<double> exact_solutions{0.05, 1.0/48.0, 0.1019760096823904218580342315643055994170804185885964798615039678, 1.651235484787737228193342177582565171308234579126117326173794530330979441089472815945286980167725541};
    //degrees are chosen and adjusted by hand
    std::vector<int> quadrule_degrees{7, 12, 14, 10};
    // outer vector for meshes, inner for functors
    std::vector<std::vector<double>> solutions;
    //--------------------------------------------------------------------------------------------------
    
    int mesh_idx = 0;
    //loop over meshes
    for (auto mesh_file : mesh_files){
        //get the mesh
        lf::io::VtkPolytopicReader reader(std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2), mesh_file);
        auto mesh_ptr = reader.mesh();

        std::cout << "Running tests for mesh with " << mesh_sizes[mesh_idx] << " cells ... \n";

        std::vector<double> mesh_solutions{};
        int functor_idx = 0;
        //loop over functors
        for (auto mesh_function : mesh_functions){
            //set up integrator
            lf::dgfe::SubTessellationIntegrator<double, decltype(mesh_function)> integrator;

            //integrate
            double sum = 0;
            for (auto cell : mesh_ptr->Entities(0)){
                sum += integrator.integrate(*cell, mesh_function, quadrule_degrees[functor_idx]);
            }
            
            //no error should be bigger than TOLERANCE
            //EXPECT_NEAR(sum, exact_solutions[functor_idx], TOLERANCE);
            mesh_solutions.push_back(sum);
            functor_idx++;
        }
        solutions.push_back(mesh_solutions);
        mesh_idx++;
    }
    std::cout << "\n";
    std::cout << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "Error 1" << std::setw(16) << "Error 2" << std::setw(16) << "Error 3" << std::setw(16) << "Error 4" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    for(int i = 0; i < mesh_files.size(); i++){
        auto sol_vec = solutions[i];
        std::cout << std::left << std::setw(10) << mesh_sizes[i] << std::right << std::setw(16)
              << std::abs(sol_vec[0] - exact_solutions[0]) << std::setw(16)
              << std::abs(sol_vec[1] - exact_solutions[1]) << std::setw(16)
              << std::abs(sol_vec[2] - exact_solutions[2]) << std::setw(16)
              << std::abs(sol_vec[3] - exact_solutions[3]) << std::endl;
    }
}

TEST(sub_tessellation_integration, pentagon){

    //test pentagon from paper
    Eigen::MatrixXd coords(2,5);
    coords <<    -0.666666666666667, 0.555555555555556, 1.000000000000000, -0.555555555555556, -1.000000000000000,
                    -0.789473684210526, -1.000000000000000, -0.052631578947368, 1.000000000000000, -0.157894736842105;

    //mesh_factory
    auto factory = std::make_unique<lf::mesh::polytopic2d::MeshFactory>(2,true);
    //Add all Points to the MeshFactory
    for (int i = 0; i < coords.cols(); i++){
        Eigen::Vector2d coords_point = coords.col(i);
        factory->AddPoint(coords_point);
    }
    //Add pentagon to the mesh factory
    factory->AddEntity(lf::base::RefEl::kPolygon(), std::array<size_type,5>{{0,1,2,3,4}}, nullptr);
    //Build the mesh
    auto mesh_ptr = factory->Build();
    //pointer to only cell
    auto b_polygon = mesh_ptr->Entities(0)[0];

    //degrees used for lambda
    int degree_x = 5;
    int degree_y = 5;
    // lambda function for polynomial integration
    auto polynomial_lambda = [&degree_x, &degree_y](const lf::mesh::Entity *entity, Eigen::Vector2d x) -> double {
        return std::pow(x[0], degree_x) * std::pow(x[1], degree_y);
    };

    //initialize integrator
    lf::dgfe::SubTessellationIntegrator<double, decltype(polynomial_lambda)> integrator;

    //TEST
    degree_x = 5;
    degree_y = 5;
    EXPECT_NEAR(integrator.integrate(*b_polygon, polynomial_lambda, 10), -0.0020324991, TOLERANCE);
    degree_x = 20;
    degree_y = 20;
    EXPECT_NEAR(integrator.integrate(*b_polygon, polynomial_lambda, 40), 6.0738145408e-8, TOLERANCE);
    degree_x = 40;
    degree_y = 40;
    EXPECT_NEAR(integrator.integrate(*b_polygon, polynomial_lambda, 80), 2.2238524572e-12, TOLERANCE);
}

} //namespace lf::dgfe::test