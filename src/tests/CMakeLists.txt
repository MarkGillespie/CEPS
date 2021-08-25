cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# https://cliutils.gitlab.io/modern-cmake/chapters/testing/googletest.html
add_subdirectory("${PROJECT_SOURCE_DIR}/deps/googletest" "deps/googletest")

# "optional"
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)
set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# This will allow you to quickly and simply add tests. Feel free to adjust to suit your
# needs. If you haven't seen it before, ARGN is "every argument after the listed ones".
# ADD YOUR LIBRARIES HERE
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main polyscope geometry-central CEPS)
    target_include_directories(${TESTNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../include" "${CMAKE_CURRENT_SOURCE_DIR}/../" "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

package_add_test(test test.cpp)
