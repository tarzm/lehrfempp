include(GoogleTest)

set(src
  integration_tests.cc
)

add_executable(lf.dgfe.test ${src})

target_link_libraries(lf.dgfe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::gtest_main lf.fe lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.polytopic2d lf.io lf.io.test_utils lf.mesh.utils lf.mesh.test_utils lf.uscalfe lf.fe.test_utils lf.dgfe)
target_compile_features(lf.dgfe.test PUBLIC cxx_std_17)
set_target_properties(lf.dgfe.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.dgfe.test)