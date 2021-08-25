#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include "Utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {

// Compute the singularities of an n-vector field
// from geometrycentral/surface/direction_fields.h

// TODO: contribute back to geometry central
VertexData<int> computeVertexIndex(ManifoldSurfaceMesh& mesh,
                                   IntrinsicGeometryInterface& geo,
                                   const FaceData<Vector2>& directionField,
                                   int nSym);

std::vector<std::pair<Vertex, double>>
lumpCones(ManifoldSurfaceMesh& mesh,
          const std::vector<std::pair<Vertex, double>>& cones);

// Doubles a mesh and geometry, gluing two copies of the mesh along their
// boundary to produce a single mesh without boundary.
// Returns the new mesh, the new geometry, the doubling map (taking each vertex
// to its twin on the other copy of the original mesh), and a list of boundary
// vertices. The doubling map sends boundary vertices to themselves
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
           std::unique_ptr<VertexPositionGeometry>, VertexData<Vertex>,
           std::vector<Vertex>>
doubleMesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);

// Copy a VertexPositionGeometry to a new ManifoldSurfaceMesh and
// EdgeLengthGeometry
std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);

// Copy an EdgeLengthGeometry to a new ManifoldSurfaceMesh and
// EdgeLengthGeometry
std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, EdgeLengthGeometry& geo);

// Count the number of flipped and zero-area triangles in a parameterization
std::pair<size_t, size_t> checkTriangleOrientations(ManifoldSurfaceMesh& mesh,
                                                    CornerData<Vector2>& uv);

namespace ImplementationDetails {

// Doubles a mesh, gluing two copies of the mesh along their
// boundary to produce a single mesh without boundary.
// Returns the new mesh, the parent vertex in the original mesh for each
// doubled vertex, and a list of boundary vertices Beware - the gluing
// involves tricky halfedge manipulations
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, VertexData<Vertex>,
           std::vector<Vertex>>
doubleMesh(ManifoldSurfaceMesh& mesh);

// Take in a list of faces, and an upper bound on the number vertices.
// Reindex the faces to remove any unused vertices and return the map
// sending an old index to its new compressed value Stolen from
// geometrycentral/surface/simple_polygon_mesh.cpp
std::vector<size_t> stripUnusedVertices(std::vector<std::vector<size_t>>& faces,
                                        size_t nV);

// Exact predicates from Jonathan Shewchuk's predicates.c
extern "C" {
void exactinit();
double orient2d(double* pa, double* pb, double* pc);
}

// Determine the orientation of triangle pqr. Returns a positive number if the
// triangle is positively oriented, a negative number if it is negatively
// oriented, and 0 if the triangle has 0 area. (Just a wrapper around Shewchuk's
// orient2d).
double orientation(const Vector2& p, const Vector2& q, const Vector2& r);


// Computes the smallest eigenvector of M^-1*A orthogonal to the all 1s vector
Vector<std::complex<double>> eig(SparseMatrix<std::complex<double>>& A,
                                 const SparseMatrix<std::complex<double>>& M,
                                 double tol);

} // namespace ImplementationDetails
} // namespace CEPS
