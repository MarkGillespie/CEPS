#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/SurfaceProjectiveParameterizationQuantity.h"
#include "polyscope/surface_mesh.h"

#include "CommonRefinement.h"
#include "ConePlacement.h"
#include "Layout.h"
#include "Optimizer.h"
#include "Tracing.h"
#include "Triangulation.h"

namespace CEPS {

// Map is stored via projective coordinates at mesh corners
struct ParameterizationResult {
    SimplePolygonMesh mesh;
    std::vector<Vector3> param;
    VertexData<int> parentMap;
    SparseMatrix<double> interpolationMatrix;
    size_t nFlippedTriangles;
    size_t nZeroAreaTriangles;
};

// Place cones using the greedy algorithm from ConePlacement.cpp
// Then uniformize, setting u=0 at cone vertices and boundary vertices.
// If viz is true, register the resulting textured mesh with polyscope.
// If checkInjectivity is true, check the intrinsically-flattened mesh for
// flipped triangles and report the results.
// uTol is the maximum allowed area distortion in the greedy cone placement
ParameterizationResult parameterizeWithGreedyCones(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, bool viz = false,
    bool checkInjectivity = false, double uTol = 5, bool verbose = false);

// Uniformize, setting u=0 at given cone vertices and any additional
// boundary vertices.
// If viz is true, register the resulting textured mesh with polyscope.
// If checkInjectivity is true, check the intrinsically-flattened mesh for
// flipped triangles and report the results.
ParameterizationResult parameterizeWithGivenCones(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
    const std::vector<Vertex>& coneVertices, bool viz = false,
    bool checkInjectivity = false, bool verbose = false);

// Uniformize, imposing the prescribed scale factors and curvatures.
// A vertex may not have both a prescribed scale factor and a prescribed
// curvature.
// Faces indicated in the imaginaryFaceIndices list will be removed from the
// final mesh (useful, e.g. for removing filled boundary loops).
// If viz is true, register the resulting textured mesh with polyscope.
// If checkInjectivity is true, check the intrinsically-flattened mesh for
// flipped triangles and report the results.
ParameterizationResult
parameterize(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
             std::vector<std::pair<Vertex, double>> prescribedScaleFactors,
             std::vector<std::pair<Vertex, double>> prescribedCurvatures,
             std::vector<size_t> imaginaryFaceIndices, bool viz = false,
             bool checkInjectivity = false, bool verbose = false);

// The core uniformization algorithm.
// The input mesh must have no boundary.
// Uniformizes the input mesh while imposing the prescribed scale factors and
// curvatures.
// All faces except for those in frontFaceIndices are deleted from the final
// mesh.
ParameterizationResult parameterizeHelper(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
    std::vector<std::pair<Vertex, double>> prescribedScaleFactors,
    std::vector<std::pair<Vertex, double>> prescribedCurvatures,
    std::set<size_t> frontFaceIndices, bool viz = false,
    bool checkInjectivity = false, bool verbose = false);

} // namespace CEPS
