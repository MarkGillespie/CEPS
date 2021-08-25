#include "PartiallyDecoratedTriangulation.h"

namespace CEPS {

PartiallyDecoratedTriangulation::PartiallyDecoratedTriangulation(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo)
    : Triangulation(mesh, geo) {
    initializeHorocycleData();
}

PartiallyDecoratedTriangulation::PartiallyDecoratedTriangulation(
    const Triangulation& otherTri, bool clearFlips)
    : Triangulation(otherTri, clearFlips) {
    initializeHorocycleData();
}

void PartiallyDecoratedTriangulation::initializeHorocycleData() {
    originalHorocyclicArcLengths = CornerData<double>(*mesh);

    for (Corner c : mesh->corners())
        originalHorocyclicArcLengths[c] = computeOriginalHorocyclicArcLength(c);

    invScaleFactors = VertexData<double>(*mesh);
    for (Vertex v : mesh->vertices())
        invScaleFactors[v] = computeInvScaleFactor(v);
}

size_t PartiallyDecoratedTriangulation::flipToDelaunay(GeometryType gType) {
    std::deque<Edge> edgesToCheck;
    EdgeData<char> inQueue(*mesh, true);
    for (Edge e : mesh->edges()) {
        edgesToCheck.push_back(e);
    }

    size_t nFlips = 0;
    while (!edgesToCheck.empty()) {

        // Get the top element from the queue of possibily non-Delaunay
        // edges
        Edge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e] = false;

        if (!isDelaunay(e)) {
            bool wasFlipped = flipEdge(e, gType);
            if (!wasFlipped) {
                std::cerr << "I should never want to do this, right?" << vendl;
                continue;
            }

            // Handle the aftermath of a flip
            nFlips++;

            // Add neighbors to queue, as they may need flipping now
            std::vector<Edge> neighboringEdges = {
                e.halfedge().next().edge(), e.halfedge().next().next().edge(),
                e.halfedge().twin().next().edge(),
                e.halfedge().twin().next().next().edge()};

            for (Edge nE : neighboringEdges) {
                if (!inQueue[nE]) {
                    edgesToCheck.emplace_back(nE);
                    inQueue[nE] = true;
                }
            }
        }
    }

    updateEdgeLengths();

    return nFlips;
}

void PartiallyDecoratedTriangulation::setScaleFactors(
    const VertexData<double>& u) {
    logScaleFactors = u;

    for (Vertex v : mesh->vertices())
        invScaleFactors[v] = computeInvScaleFactor(v);

    updateEdgeLengths();
}

void PartiallyDecoratedTriangulation::setVertexScaleFactor(Vertex v, double u) {
    logScaleFactors[v] = u;
    invScaleFactors[v] = computeInvScaleFactor(v);
    for (Edge e : v.adjacentEdges())
        geo->inputEdgeLengths[e] = currentEdgeLength(e);
    geo->refreshQuantities();
}


void PartiallyDecoratedTriangulation::updateEdgeLengths() {
    for (Edge e : mesh->edges()) {
        geo->inputEdgeLengths[e] = currentEdgeLength(e);
    }
    geo->refreshQuantities();
}

// Returns true if the vertex has finite scale factor and false otherwise
bool PartiallyDecoratedTriangulation::finite(Vertex v) const {
    return logScaleFactors[v] != std::numeric_limits<double>::infinity();
}

// Returns true if none of the vertices of the edge have infinite scale
// factors and false otherwise.
bool PartiallyDecoratedTriangulation::finite(Edge e) const {
    return finite(e.halfedge().vertex()) &&
           finite(e.halfedge().twin().vertex());
}

// Returns true if none of the vertices of the face have infinite scale
// factors and false otherwise. Assumes f is a triangle
bool PartiallyDecoratedTriangulation::finite(Face f) const {
    return finite(f.halfedge().vertex()) &&
           finite(f.halfedge().next().vertex()) &&
           finite(f.halfedge().next().next().vertex());
}

size_t PartiallyDecoratedTriangulation::degE(Vertex v) const {
    if (!finite(v)) throw_verbose_runtime_error("v = vInfinity");
    size_t deg = 0;

    // Loop counts with multiplicity
    for (Vertex w : v.adjacentVertices()) {
        if (finite(w)) deg++;
    }
    return deg;
}

size_t PartiallyDecoratedTriangulation::degF(Vertex v) const {
    if (!finite(v))
        throw_verbose_runtime_error(
            "v = vStd::Numeric_Limits<Double>::Infinity()");
    size_t deg = 0;
    for (Face f : v.adjacentFaces()) {
        if (finite(f)) deg++;
    }
    return deg;
}

bool PartiallyDecoratedTriangulation::isDelaunay(Edge e, double tol) const {
    if (incidentOnDegreeOneVertex(e)) return true;

    double oldFlipQ, newFlipQ;
    std::tie(oldFlipQ, newFlipQ) = delaunayFlipQuantity(e);

    return abs(oldFlipQ - newFlipQ) < tol || oldFlipQ >= newFlipQ;
}

bool PartiallyDecoratedTriangulation::isEssential(Edge e, double tol) const {
    double flipQ = delaunayFlipQuantity(e).first;
    return abs(flipQ) > tol;
}


void PartiallyDecoratedTriangulation::checkEdgeLengthsNaN() const {
    for (Edge e : mesh->edges()) {
        double eLen = currentEdgeLength(e);
        if (std::isnan(eLen)) {
            throw_verbose_runtime_error("Edge length NaN");
        } else if (eLen == 0) {
            throw_verbose_runtime_error("Edge length 0");
        } else if (finite(e) &&
                   eLen == std::numeric_limits<double>::infinity()) {
            throw_verbose_runtime_error(
                "Edge length std::numeric_limits<double>::infinity()");
        }
    }
}

bool PartiallyDecoratedTriangulation::flipEdge(Edge e, GeometryType gType) {
    double newLength;
    switch (gType) {
    case GeometryType::EUCLIDEAN:
        newLength = euclideanFlippedEdgeLength(e);
        break;
    case GeometryType::HYPERBOLIC:
        newLength = ptolemyFlippedEdgeLength(e);
        break;
    }

    Halfedge he = e.halfedge();

    size_t newNormalCoord = flippedNormalCoordinate(e);
    std::array<size_t, 2> newRoundabouts =
        flippedRoundabouts(he, newNormalCoord);

    bool preventSelfEdges = false;
    bool wasFlipped       = mesh->flip(e, preventSelfEdges);

    if (wasFlipped) {
        originalEdgeLengths[e] = newLength;
        normalCoordinates[e]   = newNormalCoord;

        roundaboutIndices[he]        = newRoundabouts[0];
        roundaboutIndices[he.twin()] = newRoundabouts[1];

        std::vector<Corner> neighboringCorners = {
            e.halfedge().corner(),
            e.halfedge().next().corner(),
            e.halfedge().next().next().corner(),
            e.halfedge().twin().corner(),
            e.halfedge().twin().next().corner(),
            e.halfedge().twin().next().next().corner()};

        for (Corner nC : neighboringCorners) {
            originalHorocyclicArcLengths[nC] =
                computeOriginalHorocyclicArcLength(nC);
        }

    } else {
        std::cerr << "Huh? why weren't you flipped?" << vendl;
    }

    return wasFlipped;
}

double PartiallyDecoratedTriangulation::computeInvScaleFactor(Vertex v) const {
    if (logScaleFactors[v] == std::numeric_limits<double>::infinity()) {
        return 0;
    } else {
        return exp(-logScaleFactors[v]);
    }
}

double PartiallyDecoratedTriangulation::computeOriginalHorocyclicArcLength(
    Corner c) const {

    /**
     *        / |\
     *       /  c \
     * he, e2      e1
     *     /        \
     *   |/          \
     *    ---- e0 ---->
     */
    Halfedge he = c.halfedge();
    Edge e2     = he.edge();
    Edge e0     = he.next().edge();
    Edge e1     = he.next().next().edge();

    // Gather lengths
    double e0L = originalEdgeLengths[e0];
    double e1L = originalEdgeLengths[e1];
    double e2L = originalEdgeLengths[e2];

    return e0L / (e1L * e2L);
}

VertexData<double> PartiallyDecoratedTriangulation::angleSums() const {
    geo->requireCornerAngles();
    const CornerData<double>& cornerAngle = geo->cornerAngles;
    VertexData<double> angleSum(*mesh, 0);

    for (Face f : mesh->faces()) {
        if (!finite(f)) continue;
        for (Corner c : f.adjacentCorners())
            angleSum[c.vertex()] += cornerAngle[c];
    }
    geo->unrequireCornerAngles();

    return angleSum;
}

bool PartiallyDecoratedTriangulation::satisfiesTriangleInequality() const {
    // Check triangle inequality
    for (Halfedge he : mesh->interiorHalfedges()) {

        if (!finite(he.face())) continue;

        double l1 = currentEdgeLength(he.edge());
        double l2 = currentEdgeLength(he.next().edge());
        double l3 = currentEdgeLength(he.next().next().edge());

        // TODO: fix epsilon
        if (l1 > l2 + l3 + 1e-5) {
            std::cout << "flip Quantity: "
                      << delaunayFlipQuantity(he.edge()).first << vendl;
            std::cout << "l1: " << l1 << "\t l2+l3: " << l2 + l3
                      << "\texcess: " << l1 - (l2 + l3) << "\t on edge "
                      << he.edge() << vendl;
            return false;
        }
    }

    return true;
}

SparseMatrix<double> PartiallyDecoratedTriangulation::partialCotanLaplacian() {
    SparseMatrix<double> L(mesh->nVertices(), mesh->nVertices());
    std::vector<Eigen::Triplet<double>> T;

    // if (!isHyperbolicDelaunay()) {
    //     throw_verbose_runtime_error(
    //         "I refuse to build laplacians for non delaunay meshes");
    // }

    bool careful = true; // TODO:?!
    if (careful && !satisfiesTriangleInequality()) {
        throw_verbose_runtime_error(
            "This should never happen. In any case, I refuse "
            "to build laplacians for meshes which do not "
            "satisfy the triangle inequality.");
    }

    geo->requireHalfedgeCotanWeights();
    geo->requireFaceAreas();
    HalfedgeData<double> cotanWeight = geo->halfedgeCotanWeights;
    VertexData<size_t> vertexIndex   = mesh->getVertexIndices();
    FaceData<double> area            = geo->faceAreas;

    for (Vertex v : mesh->vertices()) {
        if (!finite(v)) {
            size_t i = vertexIndex[v];
            T.emplace_back(i, i, 1);
        }
    }

    for (Edge eij : mesh->edges()) {
        double weight = 0;

        Halfedge hij = eij.halfedge();
        Halfedge hji = eij.halfedge().twin();

        // if (finite(hij.face())) weight += cotanWeight[hij];
        // if (finite(hji.face())) weight += cotanWeight[hji];

        if (finite(hij.face()) && area[hij.face()] > 0)
            weight += cotanWeight[hij];

        verbose_assert(weight == weight, "NaN cotan weight");

        double oldWeight = weight;
        if (finite(hji.face()) && area[hji.face()] > 0)
            weight += cotanWeight[hji];

        // if (weight != weight) {
        //     std::cout << "oldWeight: " << oldWeight
        //          << "\t cotanWeight: " << cotanWeight[hji] << end;
        // }

        verbose_assert(weight == weight, "NaN cotan weight");

        size_t i = vertexIndex[hij.vertex()];
        size_t j = vertexIndex[hji.vertex()];

        T.emplace_back(i, i, weight);
        T.emplace_back(j, j, weight);
        T.emplace_back(i, j, -weight);
        T.emplace_back(j, i, -weight);
    }

    geo->unrequireHalfedgeCotanWeights();
    geo->unrequireFaceAreas();

    L.setFromTriplets(T.begin(), T.end());
    return L;
}

std::pair<double, double>
PartiallyDecoratedTriangulation::delaunayFlipQuantity(Edge e) const {
    bool verbose = false;

    //======== Flip Quantity =========
    /*
     *         v2
     *        / |\
     *       / α  \
     *     e2      e1
     *     /        \
     *   |/ β      γ \
     * v3 ---- e -----> v1
     *    \ β'     γ'/|
     *     \        /
     *     e3      e4
     *       \ α' /
     *        \| /
     *         v4
     */

    Halfedge he  = e.halfedge();
    double alpha = originalHorocyclicArcLengths[he.next().next().corner()];
    double beta  = originalHorocyclicArcLengths[he.corner()];
    double gamma = originalHorocyclicArcLengths[he.next().corner()];

    double alphaPrime =
        originalHorocyclicArcLengths[he.twin().next().next().corner()];
    double betaPrime  = originalHorocyclicArcLengths[he.twin().next().corner()];
    double gammaPrime = originalHorocyclicArcLengths[he.twin().corner()];

    double s1 = invScaleFactors[he.next().vertex()];
    double s2 = invScaleFactors[he.next().next().vertex()];
    double s3 = invScaleFactors[he.vertex()];
    double s4 = invScaleFactors[he.twin().next().next().vertex()];

    if (verbose) {
        std::cout << "s1: " << s1 << "\ts2: " << s2 << "\ts3: " << s3
                  << "\ts4: " << s4 << vendl;
        std::cout << "alpha: " << alpha << "\tbeta: " << beta
                  << "\tgamma: " << gamma << vendl
                  << "\talphaPrime: " << alphaPrime
                  << "\tbetaPrime: " << betaPrime
                  << "\tgammaPrime: " << gammaPrime << vendl;
    }

    double eFlipQ = s3 * (beta + betaPrime) + s1 * (gamma + gammaPrime) -
                    s2 * alpha - s4 * alphaPrime;


    //======== Flipped Flip Quantity =========

    /*
     *         v2            |          v2
     *        / |\                     / |\
     *       /    \          |        /    \
     *     e2      e1               e2  |\  e1
     *     /        \        |      / β | β' \
     *   |/          \            |/    |     \
     * v3 ---- e -----> v1   |  v3 α    e5  α'  v1
     *    \          /|            \    |     /|
     *     \        /        |      \ γ | γ' /
     *     e3      e4               e3  |   e4
     *       \    /          |        \    /
     *        \| /                     \| /
     *         v4            |          v4
     *      original                  flipped
     */


    Edge e1 = e.halfedge().next().edge();
    Edge e2 = e.halfedge().next().next().edge();
    Edge e3 = e.halfedge().twin().next().edge();
    Edge e4 = e.halfedge().twin().next().next().edge();

    // Gather lengths
    double e0L = originalEdgeLengths[e];
    double e1L = originalEdgeLengths[e1];
    double e2L = originalEdgeLengths[e2];
    double e3L = originalEdgeLengths[e3];
    double e4L = originalEdgeLengths[e4];

    double e5L = (fmax(e2L, e4L) / e0L) * fmin(e2L, e4L) +
                 (fmax(e1L, e3L) / e0L) * fmin(e1L, e3L);

    alpha      = e5L / (e2L * e3L);
    alphaPrime = e5L / (e1L * e4L);
    beta       = e3L / (e5L * e2L);
    betaPrime  = e4L / (e5L * e1L);
    gamma      = e2L / (e3L * e5L);
    gammaPrime = e1L / (e4L * e5L);

    he = e.halfedge();
    s1 = invScaleFactors[he.next().vertex()];
    s2 = invScaleFactors[he.next().next().vertex()];
    s3 = invScaleFactors[he.vertex()];
    s4 = invScaleFactors[he.twin().next().next().vertex()];

    if (verbose) {
        std::cout << "s1: " << s1 << "\ts2: " << s2 << "\ts3: " << s3
                  << "\ts4: " << s4 << vendl;
        std::cout << "alpha: " << alpha << "\tbeta: " << beta
                  << "\tgamma: " << gamma << vendl
                  << "\talphaPrime: " << alphaPrime
                  << "\tbetaPrime: " << betaPrime
                  << "\tgammaPrime: " << gammaPrime << vendl;
    }

    double flippedFlipQ = s2 * (beta + betaPrime) + s4 * (gamma + gammaPrime) -
                          s3 * alpha - s1 * alphaPrime;

    return std::make_pair(eFlipQ, flippedFlipQ);
}

void PartiallyDecoratedTriangulation::flipToLogHyperbolicDelaunay(
    bool verbose) {


    VertexData<size_t> vIdx = mesh->getVertexIndices();
    CornerData<size_t> cIdx = mesh->getCornerIndices();

    EdgeData<double> logEdgeLength(*mesh);
    CornerData<double> horocyclicArcLength(*mesh);
    for (Edge e : mesh->edges()) {
        logEdgeLength[e] = log(originalEdgeLengths[e]);
    }

    auto lComputeHorocyclicArcLength = [&](Corner c) {
        Halfedge he = c.halfedge();
        Edge e2     = he.edge();
        Edge e0     = he.next().edge();
        Edge e1     = he.next().next().edge();

        // Gather lengths
        double logE0L = logEdgeLength[e0];
        double logE1L = logEdgeLength[e1];
        double logE2L = logEdgeLength[e2];

        // return e0L / (e1L * e2L);
        return exp(logE0L - logE1L - logE2L);
    };
    for (Corner c : mesh->corners()) {
        horocyclicArcLength[c] = lComputeHorocyclicArcLength(c);
    }

    auto lComputeFlipQuantity = [&](Edge e) {
        Halfedge he = e.halfedge();

        // clang-format off
        double alpha = horocyclicArcLength[he.next().next().corner()];
        double beta  = horocyclicArcLength[he.corner()];
        double gamma = horocyclicArcLength[he.next().corner()];

        double alphaPrime = horocyclicArcLength[he.twin().next().next().corner()];
        double betaPrime  = horocyclicArcLength[he.twin().next().corner()];
        double gammaPrime = horocyclicArcLength[he.twin().corner()];
        // clang-format on

        // return s3 * (beta + betaPrime) + s1 * (gamma + gammaPrime) - s2 *
        // alpha - s4 * alphaPrime;

        double s1 = computeInvScaleFactor(he.next().vertex());
        double s2 = computeInvScaleFactor(he.next().next().vertex());
        double s3 = computeInvScaleFactor(he.vertex());
        double s4 = computeInvScaleFactor(he.twin().next().next().vertex());
        if (s1 == 0 && s2 == 0 && s3 == 0 && s4 == 0) {
            s1 = 1;
            s2 = 1;
            s3 = 1;
            s4 = 1;
        }

        return s1 * (gamma + gammaPrime) - s2 * alpha +
               s3 * (beta + betaPrime) - s4 * alphaPrime;
    };

    auto lComputeFlippedEdgeLength = [&](Edge e) {
        Edge e1 = e.halfedge().next().edge();
        Edge e2 = e.halfedge().next().next().edge();
        Edge e3 = e.halfedge().twin().next().edge();
        Edge e4 = e.halfedge().twin().next().next().edge();

        // Gather lengths
        double logE0L = logEdgeLength[e];
        double logE1L = logEdgeLength[e1];
        double logE2L = logEdgeLength[e2];
        double logE3L = logEdgeLength[e3];
        double logE4L = logEdgeLength[e4];

        // return log [(e2L * e4L + e1L * e3L) / e0L] with the logSumExp trick
        // log(e^a + e^b) = log(e^a(1 + e^(b-a)))
        // So log(...) = logE2L + logE4L - logE0L + log(1 + exp(logE1L + logE3L
        // - logE2L - logE4L))

        double answer;
        if (logE2L + logE4L > logE1L + logE3L) {
            answer = logE2L + logE4L - logE0L +
                     log1p(exp(logE1L + logE3L - logE2L - logE4L));
        } else {
            answer = logE1L + logE3L - logE0L +
                     log1p(exp(logE2L + logE4L - logE1L - logE3L));
        }

        double e0L          = exp(logE0L);
        double e1L          = exp(logE1L);
        double e2L          = exp(logE2L);
        double e3L          = exp(logE3L);
        double e4L          = exp(logE4L);
        double simpleAnswer = (e2L * e4L + e1L * e3L) / e0L;

        // if (abs(answer - log(simpleAnswer)) > 1e-12) {
        //   std::cout << "Error: answer = " << answer << "\t log(simpleAnswer)
        //   = "
        //   << log(simpleAnswer) << vendl;
        // }


        return answer;
    };

    auto lFlipEdge = [&](Edge e) {
        double newLen = lComputeFlippedEdgeLength(e);

        size_t nNewIntersections;
        std::array<size_t, 2> newRoundabouts;
        Halfedge he = e.halfedge();
        if (!incidentOnDegreeOneVertex(e)) {
            nNewIntersections = flippedNormalCoordinate(e);
            newRoundabouts    = flippedRoundabouts(he, nNewIntersections);
        }

        // Rotates edge clockwise
        bool preventSelfEdges = false;
        bool flipped          = mesh->flip(e, preventSelfEdges);

        if (flipped) {
            logEdgeLength[e] = newLen;

            normalCoordinates[e] = nNewIntersections;

            roundaboutIndices[he]        = newRoundabouts[0];
            roundaboutIndices[he.twin()] = newRoundabouts[1];

            // fix other corner arc lengths
            std::vector<Corner> neighboringCorners = {
                e.halfedge().corner(),
                e.halfedge().next().corner(),
                e.halfedge().next().next().corner(),
                e.halfedge().twin().corner(),
                e.halfedge().twin().next().corner(),
                e.halfedge().twin().next().next().corner()};

            for (Corner nC : neighboringCorners) {
                horocyclicArcLength[nC] = lComputeHorocyclicArcLength(nC);
            }
            return true;
        } else {
            return false;
        }
    };

    auto lFlipIfNotDelaunay = [&](Edge e) {
        if (incidentOnDegreeOneVertex(e)) return false;
        double oldFlipQ = lComputeFlipQuantity(e);
        lFlipEdge(e);
        double newFlipQ = lComputeFlipQuantity(e);
        lFlipEdge(e);

        if (oldFlipQ < newFlipQ) {
            if (lFlipEdge(e)) {
                if (abs(oldFlipQ + newFlipQ) > 1e-4) {
                    std::cout << "old flipQ: " << oldFlipQ
                              << " newFlipQ: " << newFlipQ << vendl;
                    verbose_assert(
                        false, "flipping doesn't simply flip sign (err ~" +
                                   std::to_string(abs(oldFlipQ - newFlipQ)) +
                                   ")");
                }
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    };

    std::deque<Edge> edgesToCheck;
    EdgeData<char> inQueue(*mesh, true);
    for (Edge e : mesh->edges()) {
        edgesToCheck.emplace_back(e);
    }

    size_t nFlips = 0;
    while (!edgesToCheck.empty()) {

        // Get the top element from the queue of possibily non-Delaunay
        // edges
        Edge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e] = false;

        bool wasFlipped = lFlipIfNotDelaunay(e);

        if (!wasFlipped) continue;

        // Handle the aftermath of a flip
        nFlips++;

        // Add neighbors to queue, as they may need flipping now
        std::vector<Edge> neighboringEdges = {
            e.halfedge().next().edge(), e.halfedge().next().next().edge(),
            e.halfedge().twin().next().edge(),
            e.halfedge().twin().next().next().edge()};

        for (Edge nE : neighboringEdges) {
            if (!inQueue[nE]) {
                edgesToCheck.emplace_back(nE);
                inQueue[nE] = true;
            }
        }
    }
    for (Edge e : mesh->edges()) {
        originalEdgeLengths[e] = exp(logEdgeLength[e]);
    }

    for (Corner c : mesh->corners()) {
        originalHorocyclicArcLengths[c] = horocyclicArcLength[c];
    }
    updateEdgeLengths();

    if (verbose) {
        std::cout << vendl
                  << "  ... intrinsic edge flips finished. mesh is intrinsic "
                     "Delaunay after "
                  << nFlips << " flips." << vendl;
    }

    return;
}

VertexData<double>
PartiallyDecoratedTriangulation::distanceOfHorocyclesTo(Vertex v,
                                                        bool verbose) {

    // We can't check if vertices in flipTri equal v since they live in
    // different meshes. This is the next best thing.
    PartiallyDecoratedTriangulation flipTri(*this, false);

    VertexData<bool> isOriginalV(*mesh, false);
    isOriginalV[v]       = true;
    VertexData<bool> isV = isOriginalV.reinterpretTo(*flipTri.mesh);

    VertexData<double> infiniteHorocycleDistances(
        *flipTri.mesh, std::numeric_limits<double>::infinity());
    infiniteHorocycleDistances[flipTri.mesh->vertex(v.getIndex())] = 0;
    flipTri.setScaleFactors(infiniteHorocycleDistances);

    flipTri.flipToLogHyperbolicDelaunay(verbose);

    VertexData<double> distances(*flipTri.mesh);
    VertexData<double> errors(*flipTri.mesh, -2);
    VertexData<double> degree(*flipTri.mesh, 0);
    for (Vertex w : flipTri.mesh->vertices()) {
        degree[w] = w.degree();
        if (isV[w]) {
            distances[w] = 0;  // Shouldn't really be 0, but it doesn't matter.
            errors[w]    = -1; // Just to distinguish vInfinity
        } else {
            double distance = flipTri.originalEdgeLengths[w.halfedge().edge()];
            for (Halfedge he : w.incomingHalfedges()) {
                Vertex opp = he.vertex();
                // distance =
                //     fmin(distance, flipTri.originalEdgeLengths[he.edge()]);

                double distanceErr =
                    distance - flipTri.originalEdgeLengths[he.edge()];
                double relativeDistanceErr = distanceErr / distance;
                errors[w] = fmax(errors[w], abs(relativeDistanceErr));

                if (abs(relativeDistanceErr) > 1e-6) {
                    std::cout << "Distance err: " << distanceErr
                              << "\trelativeErr: " << relativeDistanceErr
                              << vendl;
                    throw_verbose_runtime_error("bad dist");
                }
            }
            distances[w] = 2 * log(distance);
        }
    }
    return distances.reinterpretTo(*mesh);
}

} // namespace CEPS
