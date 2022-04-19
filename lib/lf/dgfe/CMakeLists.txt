set(sources				
    integration.h
    integration.cc
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