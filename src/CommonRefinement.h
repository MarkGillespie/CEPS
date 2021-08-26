#pragma once

#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/utilities/utilities.h" // clamp

#include <set>

#include "DirectSum.h"
#include "Tracing.h"
#include "Triangulation.h"
#include "Utils.h"

namespace CEPS {

// Returns the common refinement, interpolated projective texture coordinates,
// and a map between vertices of Tc and vertex indices in the common refinement.
// Vertices of Ta which don't appear in the common refinement are assigned index
// -1 (This happens, e.g. when Ta has been doubled, in which case we only
// include the front vertices in the common refinement)
std::tuple<SimplePolygonMesh, std::vector<Vector3>, VertexData<int>>
computeCommonRefinement(Triangulation& Ta, Triangulation& Tb, Triangulation& Tc,
                        const VertexData<Vector3>& vertexPositions,
                        const CornerData<Vector2>& uv,
                        const std::set<Face>& frontFaces, bool verbose = false);

// Computes the common refinement, projective texture coordinates, vertex map,
// and an interpolation matrix mapping values on Ta to the common refinement
// (the interpolation matrix has dimensions |V_{common refinement}| x |V_A|)
std::tuple<SimplePolygonMesh, std::vector<Vector3>, VertexData<int>,
           SparseMatrix<double>>
computeCommonRefinementAndMatrix(Triangulation& Ta, Triangulation& Tb,
                                 Triangulation& Tc,
                                 const VertexData<Vector3>& vertexPositions,
                                 const CornerData<Vector2>& uv,
                                 const std::set<Face>& frontFaces,
                                 bool verbose = false);

namespace ImplementationDetails {

// A doublet represents a sparse vector entry, in the same way that a triplet
// represents a sparse matrix entry
// Note: this class avoids any constructors so that it is a POD type
struct Doublet {
    size_t col;
    double val;
    bool operator==(const Doublet& t) const;
    bool operator!=(const Doublet& t) const;
};

// Scalar multiplication
Doublet operator*(const double s, const Doublet& t);

// Printing
std::ostream& operator<<(std::ostream& str, const Doublet& a);

class SparseVector {
  public:
    SparseVector();
    std::vector<Doublet> doublets;

    SparseVector operator+(const SparseVector& tl) const;
    bool operator==(const SparseVector& tl) const;
};

// Scalar multiplication
SparseVector operator*(const double s, const SparseVector& tl);
SparseVector operator*(const SparseVector& tl, const double s);

// Printing
std::ostream& operator<<(std::ostream& str, const SparseVector& a);

template <typename T1, typename T2>
std::tuple<SimplePolygonMesh, std::vector<T1>, std::vector<T2>, VertexData<int>>
computeCommonRefinement(Triangulation& Ta, Triangulation& Tb, Triangulation& Tc,
                        const VertexData<T1>& initialData,
                        const CornerData<T2>& finalData,
                        const std::vector<char>& frontFaces,
                        bool verbose = false);

// ========================================================================
//                      Helpers for SliceTri
// ========================================================================

using Line    = std::array<size_t, 2>;
using Segment = std::array<size_t, 2>;

enum class LineType { A, B, Boundary };
inline std::ostream& operator<<(std::ostream& str, const LineType& lt) {
    switch (lt) {
    case LineType::A:
        str << "Type A";
        break;
    case LineType::B:
        str << "Type B";
        break;
    case LineType::Boundary:
        str << "Boundary";
        break;
    default:
        str << "Unrecognized LineType";
        break;
    }

    return str;
}

size_t next(size_t i);
size_t prev(size_t i);

// If there is a fan edge, return its index. Otherwise, return 3
// Only returns 'strict' fan edges (i.e. there is an edge to the opposite
// vertex)
size_t getFanEdge(const std::array<size_t, 3>& edgeLengths);

// If there is a fan edge, return its index. Otherwise, return 3
// Only returns 'strict' fan edges (i.e. there is an edge to the opposite
// vertex)
size_t getFanEdge(const std::array<std::vector<size_t>, 3>& edgePts);

// Important guarantee: if there is a fan edge, lines from the opposite vertex
// to that fan edge will always be oriented away from the vertex and towards the
// fan edge. Also, these lines will be in order in the output list, i.e. lines
// that end at points with earlier barycentric coordinates in the input list
// come first
std::vector<Line> getLines(const std::array<std::vector<size_t>, 3>& edgePts);

// Computes the intersection of the line a with line b
// Returns false if they don't intersect. If they do, returns true and sets
// intersection to be their intersection;
bool intersect(Line a, Line b, const std::vector<Vector2>& position, double& tA,
               double& tB);

// Computes the intersection of the line a with line b. Assumes they intersect
// and sets intersection to be their intersection;
void unsafeIntersect(Line a, Line b, const std::vector<Vector2>& position,
                     double& tA, double& tB);

void splitIntoSegments(
    const std::vector<Line>& lines,
    const std::vector<std::vector<size_t>>& lineVertices,
    const std::vector<Vector2>& vertices, const LineType& segType,
    std::vector<std::vector<size_t>>& vertexSegments,
    std::vector<std::vector<size_t>>& vertexSegmentDirections,
    std::vector<Segment>& segments, std::vector<LineType>& segmentTypes,
    bool ordered = false);


// Reads off the crossings along the boundary of face f
// Assumes that the MeshPaths always intersect halfedges positively
template <typename T>
std::tuple<std::array<std::vector<double>, 3>,
           std::array<std::vector<std::array<T, 2>>, 3>, std::vector<T>>
readEdgeCrossings(Face f, const EdgeData<MeshPath>& paths,
                  const CornerData<T>& inputData, bool projectiveInterpolate,
                  const VertexData<double>& logScaleFactors,
                  bool verbose = false);

// Lists of barycentric coordinates must start at 0, end at 1, and be sorted
//
// At fan edges, you need to set the left value of the vertex to be
// the value at the leftmost corner incident on the vertex, and the right
// value to be the value from the rightmost corner. Otherwise boundary
// interpolation fails

// Returns: the number of vertices in the refined triangle, the face-vertex
// adjacency list of the triangle, and the interpolated values of val1 and val2
// within each face
template <typename T1, typename T2>
std::tuple<size_t, std::vector<std::vector<size_t>>,
           std::vector<std::vector<T1>>, std::vector<std::vector<T2>>>
sliceTri(std::array<std::vector<double>, 3> bary1,
         std::array<std::vector<double>, 3> bary2,
         const std::array<std::vector<std::array<T1, 2>>, 3>& val1,
         const std::array<std::vector<std::array<T2, 2>>, 3>& val2,
         const std::vector<T1>& longSideCornerVal1,
         const std::vector<T2>& longSideCornerVal2, bool verbose = false);


} // namespace ImplementationDetails
} // namespace CEPS

#include "CommonRefinement.ipp"
