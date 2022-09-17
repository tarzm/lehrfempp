# Set sources
set(SOURCES
  run_convergence.h
  #fe_experiments.cc
  #manifactured_solution.cc
  )

# Define executable
#add_executable(projects.dgfem.meshes fe_experiments ${SOURCES})
add_executable(projects.dgfem.quadratic_dirichlet quadratic_dirichlet.cc ${SOURCES})
add_executable(projects.dgfem.full full.cc ${SOURCES})
add_executable(projects.dgfem.calculate_h calculate_h.cc)


# link executable with required modules
# target_link_libraries(projects.dgfem.meshes
#   PUBLIC Boost::program_options  Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.polytopic2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.fe
#   lf.uscalfe
#   )

# link executable with required modules
target_link_libraries(projects.dgfem.quadratic_dirichlet
  PUBLIC Boost::program_options  Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.polytopic2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.fe
  lf.uscalfe
  lf.dgfe
  )

# link executable with required modules
target_link_libraries(projects.dgfem.full
  PUBLIC Boost::program_options  Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.polytopic2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.fe
  lf.uscalfe
  lf.dgfe
  )

# link executable with required modules
target_link_libraries(projects.dgfem.calculate_h
  PUBLIC Boost::program_options  Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.polytopic2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.fe
  lf.uscalfe
  lf.dgfe
  )

#set_target_properties(projects.dgfem.meshes PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.quadratic_dirichlet PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.full PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.calculate_h PROPERTIES FOLDER "projects/dgfem")