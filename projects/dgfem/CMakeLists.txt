# Set sources
set(SOURCES
  )

# Define executable
#add_executable(projects.dgfem.meshes fe_experiments ${SOURCES})
add_executable(projects.dgfem.calculate_h calculate_h.cc)
add_executable(projects.dgfem.full_bvp full_bvp.cc ${SOURCES})

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

# link executable with required modules
target_link_libraries(projects.dgfem.full_bvp
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

set_target_properties(projects.dgfem.calculate_h PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.full_bvp PROPERTIES FOLDER "projects/dgfem")