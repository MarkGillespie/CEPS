#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/SurfaceProjectiveSphericalParameterizationQuantity.h"
#include "polyscope/surface_mesh.h"

#include "CommonRefinement.h"
#include "Layout.h"
#include "PartiallyDecoratedTriangulation.h"
#include "PetscWrapper.h"

namespace CEPS {

// Map is stored via projective coordinates at mesh vertices
struct SphericalUniformizationResult {
    SimplePolygonMesh mesh;
    std::vector<std::array<double, 4>> param;
    VertexData<int> parentMap;
    SparseMatrix<double> interpolationMatrix;
};

SphericalUniformizationResult
hemisphericalUniformize(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
                        bool viz = false);

SphericalUniformizationResult sphericalUniformize(ManifoldSurfaceMesh& mesh,
                                                  VertexPositionGeometry& geo,
                                                  Vertex vInfinity = Vertex(),
                                                  bool viz         = false,
                                                  bool verbose     = false);

//== The spherical uniformization energy
double sphericalObjective(PartiallyDecoratedTriangulation& tri);
Vector<double> sphericalGradient(PartiallyDecoratedTriangulation& tri);
SparseMatrix<double> sphericalHessian(PartiallyDecoratedTriangulation& tri);

VertexData<double>
sphericalGradientVertexData(PartiallyDecoratedTriangulation& tri);

namespace ImplementationDetails {
VertexData<Vector2> layOutTriangulation(PartiallyDecoratedTriangulation& tri,
                                        Vertex vInfinity);

VertexData<Vector3> mobiusCenter(ManifoldSurfaceMesh& mesh,
                                 const VertexData<Vector3>& positions,
                                 VertexData<double> area, bool verbose = false);

std::tuple<SimplePolygonMesh, std::vector<std::array<double, 4>>,
           VertexData<int>, SparseMatrix<double>>
computeCommonRefinementAndMatrix(Triangulation& Ta, Triangulation& Tb,
                                 Triangulation& Tc,
                                 const VertexData<Vector3>& initialPositions,
                                 const VertexData<Vector3>& sphericalPositions,
                                 const std::set<Face>& frontFaces,
                                 bool verbose = false);
} // namespace ImplementationDetails
} // namespace CEPS
