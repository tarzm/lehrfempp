set(sources
    polytopic2d.h
    polygon.cc
    polygon.h
    mesh.h
    mesh.cc
    mesh_factory.h
    mesh_factory.cc
)


lf_add_library(lf.mesh.polytopic2d ${sources})
target_link_libraries(lf.mesh.polytopic2d PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.hybrid2d lf.mesh.utils)
set_target_properties(lf.mesh.polytopic2d PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
    add_subdirectory(test)
endif()
