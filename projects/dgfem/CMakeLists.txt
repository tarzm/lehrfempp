# Set sources
set(SOURCES 
  fe_experiments.cc)

# Define executable
add_executable(projects.dgfem.meshes fe_experiments ${SOURCES})


# link executable with required modules
target_link_libraries(projects.dgfem.meshes
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
  )

set_target_properties(projects.dgfem.meshes PROPERTIES FOLDER "projects/dgfem")