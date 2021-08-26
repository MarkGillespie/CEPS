#include "CommonRefinement.h"

namespace CEPS {

std::tuple<SimplePolygonMesh, std::vector<Vector3>, VertexData<int>>
computeCommonRefinement(Triangulation& Ta, Triangulation& Tb, Triangulation& Tc,
                        const VertexData<Vector3>& vertexPositions,
                        const CornerData<Vector2>& uv,
                        const std::set<Face>& frontFaces, bool verbose) {

    CornerData<Vector3> homogeneousUV(*Tc.mesh);
    for (Corner c : Tc.mesh->corners()) {
        homogeneousUV[c] = Vector3{uv[c].x, uv[c].y, 1};
    }

    std::vector<char> isFrontFace(Ta.mesh->nFaces(), false);
    if (frontFaces.empty()) {
        for (size_t iF = 0; iF < isFrontFace.size(); ++iF) {
            isFrontFace[iF] = true;
        }
    } else {
        FaceData<size_t> fIdx = Ta.mesh->getFaceIndices();
        for (Face f : frontFaces) {
            isFrontFace[fIdx[f]] = true;
        }
    }


    SimplePolygonMesh commonRefinement;
    std::vector<Vector3> refinedVertexPositions;
    std::vector<Vector3> refinedHomogeneousUV;
    VertexData<int> refinedVertices;
    std::tie(commonRefinement, refinedVertexPositions, refinedHomogeneousUV,
             refinedVertices) =
        ImplementationDetails::computeCommonRefinement(
            Ta, Tb, Tc, vertexPositions, homogeneousUV, isFrontFace, verbose);

    // TODO: strip unused vertices from the common refinement

    commonRefinement.vertexCoordinates = refinedVertexPositions;

    return std::tie(commonRefinement, refinedHomogeneousUV, refinedVertices);
}

std::tuple<SimplePolygonMesh, std::vector<Vector3>, VertexData<int>,
           SparseMatrix<double>>
computeCommonRefinementAndMatrix(Triangulation& Ta, Triangulation& Tb,
                                 Triangulation& Tc,
                                 const VertexData<Vector3>& vertexPositions,
                                 const CornerData<Vector2>& uv,
                                 const std::set<Face>& frontFaces,
                                 bool verbose) {

    using ImplementationDetails::computeCommonRefinement;
    using ImplementationDetails::Doublet;
    using ImplementationDetails::SparseVector;

    CornerData<Vector3> homogeneousUV(*Tc.mesh);
    for (Corner c : Tc.mesh->corners()) {
        homogeneousUV[c] = Vector3{uv[c].x, uv[c].y, 1};
    }

    std::vector<char> isFrontFace(Ta.mesh->nFaces(), false);
    if (frontFaces.empty()) {
        for (size_t iF = 0; iF < isFrontFace.size(); ++iF) {
            isFrontFace[iF] = true;
        }
    } else {
        FaceData<size_t> fIdx = Ta.mesh->getFaceIndices();
        for (Face f : frontFaces) {
            isFrontFace[fIdx[f]] = true;
        }
    }

    // We compute the interpolation matrix from Ta to the common refinement by
    // linearly interpolating the |Va|x|Va| identity matrix.
    // More explicitly, we associate a row vector to each vertex of Ta which
    // contains a 1 at that vertex's index. Linearly interpolating these row
    // vectors to new vertices of the common refinement produces the desired
    // interpolation matrix.
    // We store these vectors as custom sparse vectors (defined in
    // CommonRefinement.h)
    VertexData<size_t> vIdxA = Ta.mesh->getVertexIndices();
    VertexData<SparseVector> trivialInterpolation(*Ta.mesh);
    for (Vertex v : Ta.mesh->vertices())
        trivialInterpolation[v].doublets.push_back(Doublet{vIdxA[v], 1});

    SimplePolygonMesh commonRefinement;
    std::vector<SparseVector> refinedInterpolation;
    std::vector<Vector3> refinedHomogeneousUV;
    VertexData<int> refinedVertices;
    std::tie(commonRefinement, refinedInterpolation, refinedHomogeneousUV,
             refinedVertices) =
        computeCommonRefinement(Ta, Tb, Tc, trivialInterpolation, homogeneousUV,
                                isFrontFace, verbose);

    // Eventually we want to strip unused vertices, since for meshes where we
    // filter out faces (i.e. doubled meshes), we are left with unused vertices
    // It's most convenient to strip them out here before we build the
    // interpolation matrix
    // We essentially just need to strip out the corresponding rows of the
    // interpolation matrix
    // So we store each vertex's index in the x coordinate of its position,
    // strip unused vertices, and then build the interpolation matrix, using
    // these x coordinates to select the correct rows
    size_t old_nV = commonRefinement.vertexCoordinates.size();
    for (size_t iV = 0; iV < old_nV; ++iV) {
        commonRefinement.vertexCoordinates[iV].x = iV;
    }
    commonRefinement.stripUnusedVertices();
    size_t nV = commonRefinement.vertexCoordinates.size(); // update nV

    // Compute vertex positions on the common refinement by interpolating the
    // original vertex positions using the interpolation matrix
    std::vector<Vector3> refinedVertexPositions;
    for (const Vector3& v : commonRefinement.vertexCoordinates) {
        size_t row = std::round(v.x);
        Vector3 pos{0, 0, 0};
        for (const Doublet& d : refinedInterpolation[row].doublets) {
            pos += d.val * vertexPositions[Ta.mesh->vertex(d.col)];
        }
        refinedVertexPositions.push_back(pos);
    }

    // Build the interpolation matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    for (size_t row = 0; row < nV; ++row) {
        size_t iV = std::round(commonRefinement.vertexCoordinates[row].x);
        for (const Doublet& d : refinedInterpolation[iV].doublets) {
            tripletList.emplace_back(row, d.col, d.val);
        }
    }
    SparseMatrix<double> interpolationMatrix(
        commonRefinement.vertexCoordinates.size(), Ta.mesh->nVertices());
    interpolationMatrix.setFromTriplets(std::begin(tripletList),
                                        std::end(tripletList));

    //== Update the vertex index map (indices have shifted since we stripped
    // unused vertices)
    // Map from common refinement -> Ta, using old indices
    // Vertices which are not present in Ta have their parent set to Vertex()
    std::vector<Vertex> parentVertex(old_nV, Vertex());
    for (Vertex v : Ta.mesh->vertices()) parentVertex[refinedVertices[v]] = v;

    // Update refinedVertices array with new indices
    for (size_t new_iV = 0; new_iV < nV; new_iV++) {
        size_t old_iV =
            std::round(commonRefinement.vertexCoordinates[new_iV].x);
        Vertex parent = parentVertex[old_iV];
        if (parent != Vertex()) {
            refinedVertices[parent] = new_iV;
        }
    }

    // Mark vertices of Ta which don't appear in the common refinement
    for (Vertex v : Ta.mesh->vertices()) {
        if (refinedVertices[v] >= static_cast<int>(nV)) {
            // invalid index - must not appear
            refinedVertices[v] = -1;
        } else {
            int new_iV = refinedVertices[v];
            int old_iV =
                std::round(commonRefinement.vertexCoordinates[new_iV].x);
            if (v != parentVertex[old_iV]) {
                refinedVertices[v] = -1;
            }
        }
    }

    commonRefinement.vertexCoordinates = refinedVertexPositions;

    return std::tie(commonRefinement, refinedHomogeneousUV, refinedVertices,
                    interpolationMatrix);
}

namespace ImplementationDetails {

// ========================================================================
//                     Sparse Vector
// ========================================================================

bool Doublet::operator==(const Doublet& t) const {
    return col == t.col && val == t.val;
}
bool Doublet::operator!=(const Doublet& t) const {
    return col != t.col || val != t.val;
}
// Scalar multiplication
Doublet operator*(const double s, const Doublet& t) {
    return Doublet{t.col, s * t.val};
}
// Printing
std::ostream& operator<<(std::ostream& str, const Doublet& a) {
    str << "< col: " << a.col << ", val: " << a.val << " >";
    return str;
}

SparseVector::SparseVector() {}

SparseVector SparseVector::operator+(const SparseVector& tl) const {
    SparseVector sum;
    sum.doublets.insert(std::end(sum.doublets), std::begin(doublets),
                        std::end(doublets));
    sum.doublets.insert(std::end(sum.doublets), std::begin(tl.doublets),
                        std::end(tl.doublets));
    return sum;
}

bool SparseVector::operator==(const SparseVector& tl) const {
    if (doublets.size() != tl.doublets.size()) return false;
    for (size_t iT = 0; iT < doublets.size(); ++iT) {
        if (doublets[iT] != tl.doublets[iT]) return false;
    }
    return true;
}

// Scalar multiplication
SparseVector operator*(const double s, const SparseVector& tl) {
    SparseVector product;
    for (const Doublet& doub : tl.doublets)
        product.doublets.emplace_back(s * doub);
    return product;
}
SparseVector operator*(const SparseVector& tl, const double s) {
    SparseVector product;
    for (const Doublet& doub : tl.doublets)
        product.doublets.emplace_back(s * doub);
    return product;
}

// Printing
std::ostream& operator<<(std::ostream& str, const SparseVector& a) {
    str << "{ ";
    for (size_t iD = 0; iD + 1 < a.doublets.size(); ++iD)
        str << a.doublets[iD] << ", ";
    if (!a.doublets.empty()) str << a.doublets[a.doublets.size() - 1];
    str << " }";
    return str;
}


// ========================================================================
//                      Helpers for SliceTri
// ========================================================================

size_t next(size_t i) { return (i + 1) % 3; }
size_t prev(size_t i) { return (i + 2) % 3; }

std::vector<Line>
getLinesTriforce(const std::array<std::vector<size_t>, 3>& edgePts) {
    std::vector<Line> lines;

    auto sideLen = [&](size_t i) { return edgePts[i].size(); };

    // Compute the lines cutting across each vertex
    for (size_t iV = 0; iV < 3; ++iV) {
        size_t nLines =
            (sideLen(next(iV)) + sideLen(prev(iV)) - sideLen(iV)) / 2;

        size_t N = sideLen(next(iV));
        for (size_t iL = 0; iL < nLines; ++iL) {
            size_t start = edgePts[next(iV)][N - iL - 1];
            size_t end   = edgePts[prev(iV)][iL];
            lines.emplace_back(Line{start, end});
        }
    }


    return lines;
}

std::vector<Line>
getLinesFan(size_t longSide,
            const std::array<std::vector<size_t>, 3>& edgePts) {

    std::vector<Line> lines;
    size_t N = edgePts[longSide].size();

    // Lines to the next edge
    for (size_t iC = 0; iC < edgePts[next(longSide)].size(); ++iC) {
        size_t start = edgePts[longSide][N - iC - 1];
        size_t end   = edgePts[next(longSide)][iC];
        lines.emplace_back(Line{start, end});
    }

    // Lines from the opposite vertex
    for (size_t iC = edgePts[prev(longSide)].size();
         iC < N - edgePts[next(longSide)].size(); ++iC) {
        // Each edge has the same index as its opposite vertex
        // So starting at longSide really means starting at the vertex opposite
        // longSide
        size_t start = longSide;
        size_t end   = edgePts[longSide][iC];
        lines.emplace_back(Line{start, end});
    }

    // Lines to the previous edge
    N = edgePts[prev(longSide)].size();
    for (size_t iC = 0; iC < edgePts[prev(longSide)].size(); ++iC) {
        size_t start = edgePts[prev(longSide)][N - iC - 1];
        size_t end   = edgePts[longSide][iC];
        lines.emplace_back(Line{start, end});
    }

    return lines;
}

// If there is a fan edge, return its index. Otherwise, return 3
size_t getFanEdge(const std::array<size_t, 3>& edgeLengths) {
    size_t longSide = 0;
    for (size_t iE = 1; iE < 3; ++iE) {
        if (edgeLengths[iE] > edgeLengths[longSide]) longSide = iE;
    }

    if (edgeLengths[longSide] >
        edgeLengths[next(longSide)] + edgeLengths[prev(longSide)]) {
        return longSide;
    } else {
        return 3;
    }
}

// If there is a fan edge, return its index. Otherwise, return 3
size_t getFanEdge(const std::array<std::vector<size_t>, 3>& edgePts) {
    std::array<size_t, 3> edgeLengths{edgePts[0].size(), edgePts[1].size(),
                                      edgePts[2].size()};
    return getFanEdge(edgeLengths);
}

std::vector<Line> getLines(const std::array<std::vector<size_t>, 3>& edgePts) {
    size_t longSide = getFanEdge(edgePts);
    if (longSide < 3) {
        return getLinesFan(longSide, edgePts);
    } else {
        return getLinesTriforce(edgePts);
    }
}

// Computes the intersection of the line a with line b
// Returns false if they don't intersect. If they do, returns true and sets
// intersection to be their intersection;
bool intersect(Line a, Line b, const std::vector<Vector2>& position, double& tA,
               double& tB) {
    double t = intersectionTime(position[a[1]], position[a[0]], position[b[1]],
                                position[b[0]]);
    if (t < 0 || t > 1) return false;

    tA = t;
    tB = intersectionTime(position[b[1]], position[b[0]], position[a[1]],
                          position[a[0]]);
    return true;
}

// Computes the intersection of the line a with line b. Assumes they intersect
// and sets intersection to be their intersection;
void unsafeIntersect(Line a, Line b, const std::vector<Vector2>& position,
                     double& tA, double& tB) {

    tA = intersectionTime(position[a[1]], position[a[0]], position[b[1]],
                          position[b[0]]);
    tB = intersectionTime(position[b[1]], position[b[0]], position[a[1]],
                          position[a[0]]);

    tA = clamp(tA, 0., 1.);
    tB = clamp(tB, 0., 1.);
}

// Segments are ordered by starting barycentric coordinate inside each line, and
// lines are ordered the same way they were in the input std::vector
// Segments in vertexSegments also have this order
void splitIntoSegments(
    const std::vector<Line>& lines,
    const std::vector<std::vector<size_t>>& lineVertices,
    const std::vector<Vector2>& vertices, const LineType& segType,
    std::vector<std::vector<size_t>>& vertexSegments,
    std::vector<std::vector<size_t>>& vertexSegmentDirections,
    std::vector<Segment>& segments, std::vector<LineType>& segmentTypes,
    bool ordered) {

    // Split up lines
    for (size_t iL = 0; iL < lines.size(); ++iL) {
        Line line                 = lines[iL];
        std::vector<size_t> verts = lineVertices[iL];
        std::vector<std::pair<double, size_t>> intersectionTimes;
        intersectionTimes.reserve(verts.size());
        for (size_t iV : verts) {
            verbose_assert(iV < vertices.size(), "vertex index too big");
            intersectionTimes.push_back(std::make_pair(
                1 - barycentricCoordinate(vertices[iV], vertices[line[0]],
                                          vertices[line[1]]),
                iV));
        }

        // std::sort sorts pairs lexicographically, so this
        // effectively sorts by barycentric coordinate
        if (!ordered) {
            std::sort(std::begin(intersectionTimes),
                      std::end(intersectionTimes));
        }

        for (size_t iV = 1; iV < intersectionTimes.size(); ++iV) {
            size_t vert     = intersectionTimes[iV].second;
            size_t prevVert = intersectionTimes[iV - 1].second;
            vertexSegments[vert].push_back(segments.size());
            vertexSegmentDirections[vert].push_back(line[0]);
            vertexSegments[prevVert].push_back(segments.size());
            vertexSegmentDirections[prevVert].push_back(line[1]);
            segments.push_back(Segment{prevVert, vert});
            segmentTypes.push_back(segType);
        }
    }
}
} // namespace ImplementationDetails
} // namespace CEPS
