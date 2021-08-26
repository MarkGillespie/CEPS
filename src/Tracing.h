#pragma once
#include "Triangulation.h"
#include "Utils.h"

#include <unordered_map>

namespace CEPS {

// Represents a point along a directed mesh edge, represented using homogeneous
// (i.e. unnormalized) barycentric coordinates
struct HalfedgePt {
    Halfedge halfedge;
    Vector2 bary;
};

// Represents a point at a mesh vertex. (Due to some Delta-complex details, it
// is useful to store a corner, rather than just a vertex, associated with the
// point).
// Ordinarily points at vertices have trivial barycentric coordinates of 1
// associated with the vertex. But since we use unnormalized barycentric
// coordinates, we actually have a nontrivial barycentric coordinate value to
// store.
struct CornerPt {
    Corner corner;
    double bary;
};

// Represent a path along a mesh. The path emanates from the corner start,
// crosses mesh edges at locations in the points list, and terminates at the
// corner end.
// The crossings in points are given in order.
// baryCoords stores the normalized barycentric coordinates of each crossing
// along the path itself (start and end always have barycentric coordinates of 0
// and 1; they are not stored in the baryCoords list)
class MeshPath {
  public:
    std::vector<HalfedgePt> points;
    std::vector<double> baryCoords;

    CornerPt start, end;
};

// Topologically trace edge e1 of T1 over T2
MeshPath traceEdge(Edge e1, Triangulation& T1, Triangulation& T2);

// Trace T1 over T2
EdgeData<MeshPath> traceTopologicalTriangulation(Triangulation& T1,
                                                 Triangulation& T2,
                                                 bool verbose = false);

// Trace T1 over T2
EdgeData<MeshPath> traceGeodesicTriangulation(Triangulation& T1,
                                              Triangulation& T2,
                                              GeometryType gType,
                                              bool verbose = false);

// Trace T1 over T2, doing as much computation on T1 as possible
// This turns out to be an easier problem in many instances
EdgeData<MeshPath> traceTransposedGeodesicTriangulation(Triangulation& T1,
                                                        Triangulation& T2,
                                                        GeometryType gType,
                                                        bool verbose = false);

// Straighten Euclidean geodesic over T
void straightenEuclidean(MeshPath& curve, Triangulation& T);

// Straighten Hyperbolic geodesic over T, according to given scale factors
// return false if straightening fails
bool straightenHyperbolic(MeshPath& curve, Triangulation& T,
                          const VertexData<double>& logScaleFactors);

// Straighten Hyperbolic geodesic over T, according to given scale factors
// return false if straightening fails
bool straightenHyperbolicIteratively(MeshPath& curve, Triangulation& T,
                                     const VertexData<double>& logScaleFactors,
                                     double tol = 1e-10);

namespace ImplementationDetails {

// Flip the MeshPath to go in the other direction
MeshPath twin(MeshPath path);

// Return an equivalent HalfedgePt where he = he.edge().halfedge()
HalfedgePt canonicalize(HalfedgePt hPt);

// Return the corresponding HalfedgePt on he.edge()'s other halfedge
HalfedgePt twin(HalfedgePt hPt);

// Take in edges of T1 traced over T2 and return edges of T2 traced over T1
EdgeData<MeshPath> transpose(const EdgeData<MeshPath>& paths, Triangulation& T1,
                             Triangulation& T2);

bool normalCoordinateTriangleInequalityViolation(Face f, Halfedge& violatingHe,
                                                 const EdgeData<size_t>& n);

// He is the twin of the halfedge in a face that the curve starts from,
// index is the index of the curve, with 1 being closest to the source
// vertex of he. The returned path starts on He, and ends at some vertex
// WARNING: 1-indexed
MeshPath traceTopologicalCurve(Halfedge he, size_t index,
                               const EdgeData<size_t>& n);

std::tuple<Halfedge, size_t> nextEdge(Halfedge he, size_t p,
                                      const EdgeData<size_t>& n);


/* Computes the Lorentz metric on R3, IE
 * v.x * w.x + v.y * w.y - v.z * w.z
 */
double lorentz(Vector3 v, Vector3 w);

/* Returns lorentz(v, v)
 */
double lorentzNorm2(Vector3 v);

/* Computes the squared distance between v and w in the Lorentz metric
 */
double lorentzDist2(Vector3 v, Vector3 w);

/* Computes the distance between v and w in the Lorentz metric
 */
double lorentzDist(Vector3 v, Vector3 w);

std::tuple<Vector3, Vector3, Vector3> computeLightConeCoords(double u, double v,
                                                             double d);

/* Takes in lightlike vectors v1, v2, v3 and distance l1, l2, l3. Computes a
 * point v4 on the light cone so that the (Minsk) distance to v1 is l1, etc
 */
Vector3 placeFourthLightConePoint(Vector3 v1, Vector3 v2, Vector3 v3, double l1,
                                  double l2, double l3);


// Returns the intersection (v.x, v.y) of the projective line a1-a2 with
// projective line b1-b2 in barycentric coordinates along line a1-a2 (i.e.
// v.x * a1 + v.y * a2 lies on line b1-b2 as well). In this version, a1, a2, b1
// and b2 must be lightlike vectors
Vector2 homogeneousProjectiveIntersectionTime(Vector3 a1, Vector3 a2,
                                              Vector3 b1, Vector3 b2);

// Returns the intersection of the projective line a1-a2 with projective line
// b1-b2 in barycentric coordinates along line a1-a2 (i.2. t*a1 + (1-t)a2 lies
// on line b1-b2 as well (up to scalar multiplication))
// a1, a2 must be lightlike vectors
double projectiveIntersectionTime(Vector3 a1, Vector3 a2, Vector3 b1,
                                  Vector3 b2);
} // namespace ImplementationDetails
} // namespace CEPS
