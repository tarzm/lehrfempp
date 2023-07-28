# Set sources
set(SOURCES
  run_convergence.h
  #fe_experiments.cc
  #manifactured_solution.cc
  )

set(DIRICHLET 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)

# Define executable
#add_executable(projects.dgfem.meshes fe_experiments ${SOURCES})
add_executable(projects.dgfem.quadratic_dirichlet quadratic_dirichlet.cc ${SOURCES})
add_executable(projects.dgfem.by_hand problem_by_hand.cc ${SOURCES})
add_executable(projects.dgfem.calculate_h calculate_h.cc)
add_executable(projects.dgfem.polytopic_from_hybrid polytopic_from_hybrid.cc ${SOURCES})
add_executable(projects.dgfem.sandbox sandbox.cc ${SOURCES})


foreach(I IN ITEMS ${DIRICHLET})
  add_executable(projects.dgfem.dirichlet_${I} dirichlet_${I}.cc ${SOURCES})
endforeach()


# link executable with required modules
foreach(I IN ITEMS ${DIRICHLET})
  target_link_libraries(projects.dgfem.dirichlet_${I}
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
endforeach()


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
target_link_libraries(projects.dgfem.by_hand
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

# link executable with required modules
target_link_libraries(projects.dgfem.polytopic_from_hybrid
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
target_link_libraries(projects.dgfem.sandbox
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

set_target_properties(projects.dgfem.quadratic_dirichlet PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.by_hand PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.calculate_h PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.polytopic_from_hybrid PROPERTIES FOLDER "projects/dgfem")
set_target_properties(projects.dgfem.sandbox PROPERTIES FOLDER "projects/dgfem")