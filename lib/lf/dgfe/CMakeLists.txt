set(sources				
    integration.h
    integration.cc
    bounding_box.h
    bounding_box.cc
    legendre_dgfe.h
    legendre_dgfe.cc
    dgfe_space.h
    dgfe_space.cc
    dgfe_providers.h
    dgfe_providers.cc
    mesh_function_dgfe.h
    mesh_function_global.h
    sub_tessellation.h
    sub_tessellation.cc
    discontinuity_penalization.h
    discontinuity_penalization.cc
)

lf_add_library(lf.dgfe ${sources})
target_link_libraries(lf.dgfe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad lf.mesh.polytopic2d lf.mesh.hybrid2d)
target_compile_features(lf.dgfe PUBLIC cxx_std_17)
set_target_properties(lf.dgfe PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()

#add_subdirectory(test_utils)