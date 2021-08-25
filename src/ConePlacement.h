#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_idt.h"

#include "Utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {

// Place cones iteratively & greedily at the vertices with the most scale
// distortion. Then run a few nonlinear Gauss-Seidel iterations to improve
// placement. Works best on a connected mesh.
std::vector<Vertex> placeCones(ManifoldSurfaceMesh& mesh,
                               IntrinsicGeometryInterface& geo, double uTol = 5,
                               size_t minNumCones           = 4,
                               size_t gaussSeidelIterations = 4);

namespace ImplementationDetails {

// Compute scale factors which give the mesh constant curvature
VertexData<double>
computeConstantCurvatureScaleFactors(ManifoldSurfaceMesh& mesh,
                                     IntrinsicGeometryInterface& geo,
                                     const SparseMatrix<double>& L);

// Compute scale factors which flatten the mesh except at the given cone
// vertices If vSkip1 or vSkip2 are set, we flatten the mesh at those cones as
// well.
VertexData<double> computeScaleFactors(ManifoldSurfaceMesh& mesh,
                                       IntrinsicGeometryInterface& geo,
                                       const SparseMatrix<double>& L,
                                       const std::vector<Vertex> cones,
                                       int vSkip1 = -1, int vSkip2 = -1);


SparseMatrix<double> intrinsicDelaunayLaplacian(ManifoldSurfaceMesh& mesh,
                                                IntrinsicGeometryInterface& geo,
                                                double mollify = 1e-6);
} // namespace ImplementationDetails
} // namespace CEPS
