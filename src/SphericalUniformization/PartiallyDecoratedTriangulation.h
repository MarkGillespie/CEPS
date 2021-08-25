#pragma once
#include <deque>

#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/elementary_geometry.h" // placeTriangleVertex

#include "FlipFormulas.h"
#include "GeometryUtils.h"
#include "Triangulation.h"
#include "Utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {

class PartiallyDecoratedTriangulation : public Triangulation {
  public:
    PartiallyDecoratedTriangulation(ManifoldSurfaceMesh& mesh,
                                    VertexPositionGeometry& geo);

    PartiallyDecoratedTriangulation(const Triangulation& otherTri,
                                    bool clearFlips);

    //== Core data
    VertexData<double> invScaleFactors;
    CornerData<double> originalHorocyclicArcLengths;

    //== Mutators

    void initializeHorocycleData();

    // Flip to intrinsic Delaunay if gType is EUCLIDEAN and ideal Delaunay if
    // gType is HYPERBOLIC
    size_t flipToDelaunay(GeometryType gType);

    // Flip to ideal Delaunay, storing geometry via log lengths
    // (Generally you can just use flipToDelaunay. This should only be important
    // when computing single-source geodesic distance)
    void flipToLogHyperbolicDelaunay(bool verbose);

    void setScaleFactors(const VertexData<double>& u);
    void setVertexScaleFactor(Vertex v, double u);

    void updateEdgeLengths(); // TODO: also update invscalefactors

    // Updates mesh, normal coordinates, roundabouts, edge lengths, and
    // horocyclic arc lengths.
    // Uses Euclidean formula if gType is EUCLIDEAN and
    // Ptolemy formula if gType is HYPERBOLIC Returns true if edge was flipped,
    // false if edge was unflippable
    bool flipEdge(Edge e, GeometryType gType);

    //== Data Accesses
    /* double currentEdgeLength(Edge e) const; */

    // Return the (non-log) inverse scale factor associated with vertex v
    double computeInvScaleFactor(Vertex v) const;

    // Returns the length of the horocyclic arc relative to the original
    // horocycles
    double computeOriginalHorocyclicArcLength(Corner c) const;

    VertexData<double> angleSums() const;

    // Returns true if the vertex has finite scale factor and false otherwise
    bool finite(Vertex v) const;

    // Returns true if none of the vertices of the edge have infinite scale
    // factors and false otherwise.
    bool finite(Edge f) const;

    // Returns true if none of the vertices of the face have infinite scale
    // factors and false otherwise. Assumes f is a triangle
    bool finite(Face f) const;

    // Numer of finite edges incident on vertex v
    size_t degE(Vertex v) const;

    // Numer of finite faces incident on vertex v
    size_t degF(Vertex v) const;

    bool isDelaunay(Edge e, double tol = 1e-8) const;

    // Essential edges of a Delaunay triangulation may not be flipped.
    // Nonessential edges may
    bool isEssential(Edge e, double tol = 1e-8) const;

    void checkEdgeLengthsNaN() const;

    bool satisfiesTriangleInequality() const;

    // Cotan laplacian of finite part of mesh
    SparseMatrix<double> partialCotanLaplacian();

    // Distance from v to v is given as 0 even though that isn't the true answer
    // because we never need it and I don't know how to compute it
    VertexData<double> distanceOfHorocyclesTo(Vertex v, bool verbose = false);

    // Delaunay flip quantity for edge e, and flip(e)
    std::pair<double, double> delaunayFlipQuantity(Edge e) const;

    void initializeRoundabouts();
};
} // namespace CEPS
