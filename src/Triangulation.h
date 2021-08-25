#pragma once
#include <deque>

#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/elementary_geometry.h" // placeTriangleVertex

#include "FlipFormulas.h"
#include "GeometryUtils.h"
#include "Utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {

class Triangulation {
  public:
    Triangulation(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);
    Triangulation(const Triangulation& otherTri, bool clearFlips = false);

    // connectivity
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    // geometry of rescaled edge lengths
    std::unique_ptr<EdgeLengthGeometry> geo;

    // Flip to intrinsic Delaunay if gType is EUCLIDEAN and ideal Delaunay if
    // gType is HYPERBOLIC
    size_t flipToDelaunay(GeometryType gType);

    void setScaleFactors(const VertexData<double>& u);

    VertexData<double> logScaleFactors;
    EdgeData<double> originalEdgeLengths;
    EdgeData<size_t> normalCoordinates;
    HalfedgeData<size_t> roundaboutIndices;
    VertexData<size_t> roundaboutDegrees;

    double currentEdgeLength(Edge e) const;
    void updateEdgeLengths();

    // Updates mesh, normal coordinates, roundabouts, and edge lengths
    // Uses Euclidean formula if gType is EUCLIDEAN and Ptolemy formula if gType
    // is HYPERBOLIC
    // Returns true if edge was flipped, false if edge was unflippable
    bool flipEdge(Edge e, GeometryType gType);

    bool incidentOnDegreeOneVertex(Edge e) const;

    bool isDelaunay(Edge e) const;

    // Delaunay flip quantity (eq 9) for edge e, and flip(e)
    std::pair<double, double> delaunayFlipQuantity(Edge e) const;

    double euclideanFlippedEdgeLength(Edge e) const;
    double ptolemyFlippedEdgeLength(Edge e) const;
    size_t flippedNormalCoordinate(Edge e) const;
    std::array<size_t, 2> flippedRoundabouts(Halfedge he,
                                             size_t flippedNormalCoord) const;

    // Helper functions to gather the data needed for all the miscellaneous flip
    // formulas
    std::array<double, 5> currentNeighborhoodEdgeLengths(Halfedge he) const;
    std::array<double, 5> currentNeighborhoodEdgeLengths(Edge e) const;
    std::array<double, 5> originalNeighborhoodEdgeLengths(Halfedge he) const;
    std::array<double, 5> originalNeighborhoodEdgeLengths(Edge e) const;
    std::array<size_t, 5> neighborhoodNormalCoordinates(Halfedge he) const;
    std::array<size_t, 5> neighborhoodNormalCoordinates(Edge e) const;

    void initializeRoundabouts();
};
} // namespace CEPS
