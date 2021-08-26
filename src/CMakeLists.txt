cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()


add_library(
  CEPS
  Triangulation.cpp
  Optimizer.cpp
  Tracing.cpp
  FlipFormulas.cpp
  Layout.cpp
  CommonRefinement.cpp
  Utils.cpp
  GeometryUtils.cpp
  Cutter.cpp
  ConePlacement.cpp
  ConeFlattening.cpp
  IO.cpp
  Clausen.cpp
  Logger.cpp
  polyscope/SurfaceProjectiveParameterizationQuantity.cpp
  ../deps/predicates/predicates.c
)
target_link_libraries(CEPS polyscope geometry-central)

add_executable(
  parameterize
  Parameterize.cpp
)
target_include_directories(parameterize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
# borrow args.hxx directly from polyscope
target_include_directories(parameterize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/polyscope/deps/args/")
target_link_libraries(parameterize polyscope geometry-central CEPS)

option(SPHERICAL_UNIFORMIZATION "Build spherical uniformizer" OFF)
option(PACKAGE_TESTS "Build the cone flattening tests" OFF)

message("Spherical uniformization: ${SPHERICAL_UNIFORMIZATION}")
message("Tests: ${PACKAGE_TESTS}")

if(SPHERICAL_UNIFORMIZATION)
  # PETSC_DIR=/home/mark/Documents/Code/petsc PETSC_ARCH=arch-linux2-c-opt all
  add_subdirectory(SphericalUniformization)
endif()

if(PACKAGE_TESTS)
  enable_testing()
  add_subdirectory(
    tests
  )
endif()
