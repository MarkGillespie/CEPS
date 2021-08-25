include("findPETSc.cmake")

add_library(
  CEPS_SPHERICAL
  SphericalUniformization.cpp
  PartiallyDecoratedTriangulation.cpp
  PetscWrapper.cpp
  polyscope/SurfaceProjectiveSphericalParameterizationQuantity.cpp
  ../../deps/polyscope/deps/stb/stb_impl.cpp
  polyscope/bindata_map_equirectangular.cpp
)
target_include_directories(CEPS_SPHERICAL PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${PETSC_INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_include_directories(CEPS_SPHERICAL PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/stb/")
target_link_libraries(CEPS_SPHERICAL polyscope geometry-central CEPS ${PETSC_LIBRARIES})

add_executable(
  spherical_uniformize
  SphericalDemo.cpp
  )
target_include_directories(spherical_uniformize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${PETSC_INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/..")
# borrow args.hxx directly from polyscope
target_include_directories(spherical_uniformize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args/")
target_link_libraries(spherical_uniformize polyscope geometry-central CEPS CEPS_SPHERICAL ${PETSC_LIBRARIES})
