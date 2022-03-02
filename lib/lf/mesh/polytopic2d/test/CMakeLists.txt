include(GoogleTest)

set(sources
    polygon_test.cc
    polygon_test.h)

add_executable(lf.mesh.polytopic2d.test ${sources})
target_link_libraries(lf.mesh.polytopic2d.test PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main
  lf.io
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.polytopic2d)
set_target_properties(lf.mesh.polytopic2d.test PROPERTIES FOLDER "lib")
target_compile_features(lf.mesh.polytopic2d.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.mesh.polytopic2d.test)