/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/dgfe/dgfe.h>

#include "run_convergence.h"

#include <fstream>
#include <iomanip>
#include <string>

void write_error_file(std::string run_name, double c_inv, int c_sigma, int num_cells, std::string error_type, double error){
    //error file
    auto c_inv_str = std::to_string(c_inv);
    c_inv_str.resize(4);
    auto c_sigma_str = std::to_string(c_sigma);
    std::string out_file_name = "measurements/" + run_name + "/" + std::to_string(num_cells) + "_" + c_inv_str
                                 + "_" + c_sigma_str + "_" + error_type + ".txt";
    std::cout << "Trying to write to " << out_file_name << "\n";
    std::ofstream out_file(out_file_name);
    out_file << error;
    out_file.close();
}

int main(int argc, char* argv[]) {

std::string run_name = argv[1];
std::cout << "\nRUN NAME " << run_name << "\n\n";



// abbreviations for types
using size_type = lf::base::size_type;
using glb_idx_t = lf::assemble::glb_idx_t;
using coord_t = Eigen::Vector2d;

// ======================================================================
  // Prep stage: provide all coefficient functions mainly through lambda
  //             functions and derived MeshFunctions.
  // ======================================================================

  // Coefficients:

  // 2x2 diffusion tensor A(x)
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return (Eigen::Matrix<double, 2, 2>() << 1.0, 0.0, 0.0 , 1.0)
        .finished();
  };
  // Wrap diffusion coefficient into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};

  // Scalar valued reaction coefficient c
  auto gamma = [](Eigen::Vector2d x) -> double {
    return (0.0);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_gamma{gamma};

  // vector valued advection coefficient
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d { return (Eigen::Matrix<double, 2, 1>() << 0.0, 0.0).finished(); };
  lf::mesh::utils::MeshFunctionGlobal mf_beta{beta};

  /* SAM_LISTING_BEGIN_1 */
  // Exact solution u
  auto u = [](Eigen::Vector2d x) -> double {
    return 1 + x[0] + x[1]*x[1];
  };
  // Has to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};

  // Gradient of exact solution
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    double den = x[0] * x[0] + x[1] + 1.0;
    return ((Eigen::Vector2d() << 1.0, 1.0).finished());
  };
  // Convert into mesh function to use for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  /* SAM_LISTING_END_1 */

  // Right-hand side source function f
  auto f = [&gamma, &u](Eigen::Vector2d x) -> double {
    return -2.0;
  };
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Dirichlet data borrowed from the known exact solution.
  auto g = [&u](const Eigen::Vector2d& x) -> double { return u(x); };
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};

  // Predicates for selecting edges on the Dirichlet boundary
  std::function<bool(const Eigen::Vector2d&)> dir_sel =
      [](const Eigen::Vector2d& x) -> bool { return true; };
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(const Eigen::Vector2d&)>>
      edge_sel_dir{dir_sel};


  // ======================================================================
  // Stage I: Definition of computational domain through coarsest mesh
  // Since this example relies on a manufactured solution tied to a particular
  // domain, using a hard-wired mesh is justified. Another example will address
  // solving a boundary value problem on a mesh read from file.
  // ======================================================================
  // The following code also illustrates the role of a MeshFactory
  // Create helper object: mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Generate nodes of the mesh
  // clang-format off
  std::array<std::array<double, 2>, 8> node_coord{
  std::array<double, 2>({0   , 0    }),
  std::array<double, 2>({1   , 0    }),
  std::array<double, 2>({0.2 , 1.0  }),
  std::array<double, 2>({0.5 , 0.0  }),
  std::array<double, 2>({0.6 , 0.5  }),
  std::array<double, 2>({0.1 , 0.5  }),
  std::array<double, 2>({0.15, 0.75 }),
  std::array<double, 2>({0.4 , 0.75 })};
  // clang-format on
  // Add nodes to the mesh via the MeshFactory object
  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(coord_t({node[0], node[1]}));
  }
  // Add plain triangles to the mesh, defined by their vertex nodes.
  // Since no particular geometry is specified, the triangles are assumed to
  // have straight edges.
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({3, 1, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({7, 6, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Create a general quadrilateral with straight edges.
  std::array<size_type, 4> quad_nodes{5, 4, 7, 6};
  Eigen::Matrix<double, 2, 4> quad_coord;
  for (int n_pt = 0; n_pt < 4; ++n_pt) {
    quad_coord(0, n_pt) = node_coord[quad_nodes[n_pt]][0];
    quad_coord(1, n_pt) = node_coord[quad_nodes[n_pt]][1];
  }
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(), std::vector<size_type>({5, 4, 7, 6}),
      std::make_unique<lf::geometry::QuadO1>(quad_coord));
  // Create parallelogram.
  std::array<size_type, 4> parg_nodes{0, 3, 4, 5};
  for (int n_pt = 0; n_pt < 4; ++n_pt) {
    quad_coord(0, n_pt) = node_coord[parg_nodes[n_pt]][0];
    quad_coord(1, n_pt) = node_coord[parg_nodes[n_pt]][1];
  }
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(), std::vector<size_type>({0, 3, 4, 5}),
      std::make_unique<lf::geometry::Parallelogram>(quad_coord));

  // Get a pointer to the coarsest mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // Print information about the coarsest mesh
  std::cout << "\t Coarsest mesh for demonstration run\n";
  //lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

  // ======================================================================
  // Stage II: Ask LehrFEM++ to create a hierarchy of nested meshes
  // ======================================================================


// Obtain a pointer to a hierarchy of nested meshes
const int reflevels = 4;
std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
    lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                            reflevels);
lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
// Ouput information about hierarchy of nested meshes
std::cout << "\t Sequence of nested meshes used in demo code\n";
multi_mesh.PrintInfo(std::cout);
// Number of levels
size_type L = multi_mesh.NumLevels();

// Vector for keeping error norms
std::vector<std::tuple<size_type, double, double>> errs{};
// Vector for keeping polytopic error norms
std::vector<std::tuple<size_type, double, double>> errs_poly{};

//define parameters
double c_inv = 0.5;
double c_sigma = 40.0;


// ############################LEVEL LOOP: Do computations on all levels
for (size_type level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Reference to current mesh
    const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

    // Preprocessing: count number of edges with different boundary conditions
    size_type no_Dirichlet_edges = 0;
    size_type no_Neumann_edges = 0;
    size_type no_impedance_edges = 0;
    // Obtain an array of boolean flags for the edges of the mesh, 'true'
    // indicates that the edge lies on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
    // Traverse the edges of the mesh and check their boundary flags and the
    // type of boundary condition
    for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
        if (bd_flags(*edge)) {
            if (edge_sel_dir(*edge)) {
                no_Dirichlet_edges++;
            }
        }
    }
    // Dimension of finite element space`
    const size_type N_dofs(dofh.NumDofs());

    // Verbose output
    std::cout << "Computations on level " << level
                << " (#Dir_ed =  " << no_Dirichlet_edges
                << ", #Neu_ed = " << no_Neumann_edges
                << ", #imp_ed = " << no_impedance_edges << "), #dof = " << N_dofs
                << std::endl;

    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

    // ----------------------------------------------------------------------
    // III: Assemble finite element Galerkin matrix
    // First the volume part for the bilinear form
    // Initialize object taking care of local computations. No selection of a
    // subset of cells is specified in this demonstration: assembly will cover
    // all cells.
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<
        double, decltype(mf_alpha), decltype(mf_gamma)>
        elmat_builder(fe_space, mf_alpha, mf_gamma);
    // Invoke assembly on cells (co-dimension = 0 as first argument)
    // Information about the mesh and the local-to-global map is passed through
    // a Dofhandler object, argument 'dofh'. This function call adds triplets to
    // the internal COO-format representation of the sparse matrix A.
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

    // ----------------------------------------------------------------------
    // IV: Right-hand side vector; has to be set to zero initially
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
    phi.setZero();
    // Assemble volume part of right-hand side vector depending on the source
    // function f.
    // Initialize object taking care of local computations on all cells.
    lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
        elvec_builder(fe_space, mf_f);
    // Invoke assembly on cells (codim == 0)
    AssembleVectorLocally(0, dofh, elvec_builder, phi);

    // ----------------------------------------------------------------------
    // III: Fixing solution components according to essential (Dirichlet)
    // boundary conditions
    if (no_Dirichlet_edges > 0) {
        // Obtain specification for shape functions on edges
        const auto* rsf_edge_p =
            fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
        LF_ASSERT_MSG(rsf_edge_p != nullptr,
                    "FE specification for edges missing");

        // Fetch flags and values for degrees of freedom located on Dirichlet
        // edges.
        auto ess_bdc_flags_values{lf::fe::InitEssentialConditionFromFunction(
            *fe_space,
            [&edge_sel_dir, &bd_flags](const lf::mesh::Entity& edge) -> bool {
            return (bd_flags(edge) && edge_sel_dir(edge));
            },
            mf_g)};
        // Eliminate Dirichlet dofs from linear system
        lf::assemble::FixFlaggedSolutionComponents<double>(
            [&ess_bdc_flags_values](glb_idx_t gdof_idx) {
            return ess_bdc_flags_values[gdof_idx];
            },
            A, phi);
    }

    /* SAM_LISTING_BEGIN_2 */
    // Assembly completed: Convert COO matrix A into CRS format using Eigen's
    // internal conversion routines.
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // Solve linear system using Eigen's sparse direct elimination
    // Examine return status of solver in case the matrix is singular
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(phi);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    // Postprocessing: Compute error norms
    // create mesh functions representing solution / gradient of solution
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            mesh, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        mesh, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    /* SAM_LISTING_END_2 */

    //write error
    auto n_cells = mesh_p->Entities(0).size();
    errs.emplace_back(n_cells, L2err, H1serr);
    write_error_to_file(n_cells, L2err, run_name);

    //############################ END HYBRID

    //############################ START POLYTOPIC

    //mesh
    auto poly_mesh_ptr = lf::mesh::polytopic2d::PolytopicFromHybrid2D(mesh_p);

    //write mesh for python drawing
    //lf::io::writeMatplotlib(*poly_mesh_ptr, "./csvs/" + std::to_string(poly_mesh_ptr->NumEntities(0)) + ".csv");

    //dgfe space
    lf::dgfe::DGFESpace dgfe_space(poly_mesh_ptr, 2);
    auto dgfe_space_ptr = std::make_shared<lf::dgfe::DGFESpace>(dgfe_space);

    // setup of mesh functions
    lf::dgfe::MeshFunctionGlobalDGFE m_a_coeff{alpha};
    lf::dgfe::MeshFunctionGlobalDGFE m_b_coeff{beta};
    lf::dgfe::MeshFunctionGlobalDGFE m_c_coeff{gamma};
    lf::dgfe::MeshFunctionGlobalDGFE m_f{f};
    lf::dgfe::MeshFunctionGlobalDGFE m_gD{f};
    lf::dgfe::MeshFunctionGlobalDGFE m_gN{f};


    // Setup l2 projection of sqrt(A) * nabla(basis)
    auto l2_projection = lf::dgfe::L2ProjectionSqrtANablaBasis<double>(dgfe_space_ptr, m_a_coeff, 20);

    // double c_inv_opt = 100.0;
    // double c_sigma_opt = 100.0;
    // double smallest_err = 100.0;

    // for (double c_inv = 0.05; c_inv < 1; c_inv += 0.1){
    //     for (double c_sigma = 10.0; c_sigma < 200.0 ; c_sigma += 10.0){
    //         //run it
    //         auto error_l2_poly = run_convergence(c_inv, c_sigma, 10, run_name, dgfe_space_ptr, l2_projection, m_a_coeff, m_b_coeff, m_c_coeff, m_gD, m_gN, m_f, m_gD);
    //         if (error_l2_poly < smallest_err){
    //             c_inv_opt = c_inv;
    //             c_sigma_opt = c_sigma;
    //             smallest_err = error_l2_poly;
    //         }
    //     }
    // }

    

    auto error_l2_poly = run_convergence(c_inv, c_sigma, 5, run_name, dgfe_space_ptr, l2_projection, m_a_coeff, m_b_coeff, m_c_coeff, m_gD, m_gN, m_f, m_gD);

    //std::cout << "Smallest error for c_inv and c_sigma: " << smallest_err << " for " << c_inv_opt << " and " << c_sigma_opt << "\n\n";

    //run it
    // auto error_l2_poly = run_convergence(c_inv, c_sigma, 10, run_name, dgfe_space_ptr, l2_projection, m_a_coeff, m_b_coeff, m_c_coeff, m_gD, m_gN, m_f, m_gD);

    errs_poly.emplace_back(n_cells, error_l2_poly, error_l2_poly);

    auto num_cells = poly_mesh_ptr->NumEntities(0);
    write_error_file(run_name, c_inv, c_sigma, num_cells, "L2_poly", error_l2_poly);
    write_error_file(run_name, c_inv, c_sigma, num_cells, "L2_hybrid", L2err);



}
// ############################ END LEVEL LOOP: Do computations on all levels


//###########HYBRID
// Output table of errors to file and terminal
std::cout << "\nHYBRID\n";
std::ofstream out_file("errors.txt");
std::cout << std::left << std::setw(10) << "N cells" << std::right << std::setw(16)
        << "L2 error" << std::setw(16) << "H1 error" << std::endl;
std::cout << "---------------------------------------------" << std::endl;
for (const auto& err : errs) {
auto [N, l2err, h1serr] = err;
out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
            << l2err << std::setw(16) << h1serr << std::endl;
std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
            << l2err << std::setw(16) << h1serr << std::endl;
}

//############POLYTOPIC
std::cout << "\n\nPOLYTOPIC with C_sigma = " << c_sigma << " and C_INV = " << c_inv << "\n" ;
std::cout << std::left << std::setw(10) << "N cells" << std::right << std::setw(16)
        << "L2 error" << std::setw(16) << "H1 error" << std::endl;
std::cout << "---------------------------------------------" << std::endl;
for (const auto& err : errs_poly) {
auto [N, l2err, h1serr] = err;
out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
            << l2err << std::setw(16) << h1serr << std::endl;
std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
            << l2err << std::setw(16) << h1serr << std::endl;
}



return 0;
}
