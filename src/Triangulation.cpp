#include "Triangulation.h"

namespace CEPS {
Triangulation::Triangulation(ManifoldSurfaceMesh& mesh_,
                             VertexPositionGeometry& geo_) {
    std::tie(mesh, geo) = copyGeometry(mesh_, geo_);

    mollifyIntrinsic(*mesh, geo->inputEdgeLengths, 1e-12);

    originalEdgeLengths = geo->inputEdgeLengths;
    logScaleFactors     = VertexData<double>(*mesh, 0);
    normalCoordinates   = EdgeData<size_t>(*mesh, 0);
    initializeRoundabouts();
}

Triangulation::Triangulation(const Triangulation& otherTri, bool clearFlips) {
    std::tie(mesh, geo) = copyGeometry(*otherTri.mesh, *otherTri.geo);

    originalEdgeLengths = otherTri.originalEdgeLengths.reinterpretTo(*mesh);
    logScaleFactors     = otherTri.logScaleFactors.reinterpretTo(*mesh);
    if (clearFlips) {
        normalCoordinates = EdgeData<size_t>(*mesh, 0);
        initializeRoundabouts();
    } else {
        normalCoordinates = otherTri.normalCoordinates.reinterpretTo(*mesh);
        roundaboutIndices = otherTri.roundaboutIndices.reinterpretTo(*mesh);
        roundaboutDegrees = otherTri.roundaboutDegrees.reinterpretTo(*mesh);
    }
}

size_t Triangulation::flipToDelaunay(GeometryType gType) {
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

                std::cout << e << vendl;
                std::cout << e.halfedge().face() << "\t"
                          << e.halfedge().twin().face() << vendl;

                geo->requireCornerAngles();
                std::cout << "corner angles: " << vendl;
                for (Corner c : e.halfedge().face().adjacentCorners())
                    std::cout << "\t" << geo->cornerAngles[c] << vendl;
                geo->unrequireCornerAngles();
                std::cout << "Edge lengths: " << vendl;
                for (Halfedge he : e.halfedge().face().adjacentHalfedges())
                    std::cout << "\t" << geo->inputEdgeLengths[he.edge()]
                              << vendl;

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

void Triangulation::setScaleFactors(const VertexData<double>& u) {
    logScaleFactors = u;
    updateEdgeLengths();
}

double Triangulation::currentEdgeLength(Edge e) const {
    double u1 = logScaleFactors[src(e)];
    double u2 = logScaleFactors[dst(e)];
    return exp((u1 + u2) / 2.) * originalEdgeLengths[e];
}

void Triangulation::updateEdgeLengths() {
    for (Edge e : mesh->edges()) {
        geo->inputEdgeLengths[e] = currentEdgeLength(e);
    }
    geo->refreshQuantities();
}

bool Triangulation::incidentOnDegreeOneVertex(Edge e) const {
    // Get half edges of first face
    Halfedge ha1 = e.halfedge();
    Halfedge ha2 = ha1.next();

    // Get halfedges of second face
    Halfedge hb1 = ha1.twin();
    Halfedge hb2 = hb1.next();

    // Check whether incident on degree 1 vertex
    return ha2 == hb1 || hb2 == ha1;
}

bool Triangulation::isDelaunay(Edge e) const {
    if (incidentOnDegreeOneVertex(e)) return true;

    double q, flippedQ;
    std::tie(q, flippedQ) = delaunayFlipQuantity(e);
    return q >= flippedQ;
}

// Delaunay flip quantity (eq 9) for edge e, and flip(e)
std::pair<double, double> Triangulation::delaunayFlipQuantity(Edge e) const {
    return delaunayCondition(currentNeighborhoodEdgeLengths(e));
}

// Rotates edge clockwise
bool Triangulation::flipEdge(Edge e, GeometryType gType) {
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
    } else {
        std::cerr << "Huh? why weren't you flipped?" << vendl;
    }

    return wasFlipped;
}

size_t Triangulation::flippedNormalCoordinate(Edge e) const {
    return flipNormalCoordinate(neighborhoodNormalCoordinates(e));
}

std::array<size_t, 2>
Triangulation::flippedRoundabouts(Halfedge he,
                                  size_t flippedNormalCoord) const {
    // Gather data needed for flip
    size_t nlk = flippedNormalCoord;

    size_t dk = roundaboutDegrees[he.next().next().vertex()];
    size_t dl = roundaboutDegrees[he.twin().next().next().vertex()];

    size_t rki = roundaboutIndices[he.next().next()];
    size_t rlj = roundaboutIndices[he.twin().next().next()];

    return flipRoundabout(neighborhoodNormalCoordinates(he), nlk, dk, dl, rki,
                          rlj);
}

double Triangulation::euclideanFlippedEdgeLength(Edge e) const {
    double currentLen = flipEuclideanLength(currentNeighborhoodEdgeLengths(e));

    double u1 = logScaleFactors[e.halfedge().next().next().vertex()];
    double u2 = logScaleFactors[e.halfedge().twin().next().next().vertex()];

    double originalLen = currentLen / exp(0.5 * (u1 + u2));

    return originalLen;
}

double Triangulation::ptolemyFlippedEdgeLength(Edge e) const {
    return flipHyperbolicLength(originalNeighborhoodEdgeLengths(e));
}

void Triangulation::initializeRoundabouts() {
    roundaboutIndices = HalfedgeData<size_t>(*mesh);
    roundaboutDegrees = VertexData<size_t>(*mesh);
    for (Vertex v : mesh->vertices()) {
        size_t iHe           = 0;
        size_t D             = v.degree();
        roundaboutDegrees[v] = D;

        // Explicitly loop over halfedges in counterclockwise order
        Halfedge he = v.halfedge();
        do {
            roundaboutIndices[he] = iHe;

            iHe = (iHe + 1) % D;
            he  = he.next().next().twin();
        } while (he != v.halfedge());
    }
}

// Helper functions to gather the data needed for all the miscellaneous flip
// formulas
std::array<double, 5>
Triangulation::currentNeighborhoodEdgeLengths(Halfedge he) const {
    // Gather current lengths
    double lij = currentEdgeLength(he.edge());
    double ljk = currentEdgeLength(he.next().edge());
    double lki = currentEdgeLength(he.next().next().edge());
    double lil = currentEdgeLength(he.twin().next().edge());
    double llj = currentEdgeLength(he.twin().next().next().edge());
    return {lij, ljk, lki, lil, llj};
}

std::array<double, 5>
Triangulation::currentNeighborhoodEdgeLengths(Edge e) const {
    return currentNeighborhoodEdgeLengths(e.halfedge());
}

std::array<double, 5>
Triangulation::originalNeighborhoodEdgeLengths(Halfedge he) const {
    // Gather original lengths
    double lij = originalEdgeLengths[he.edge()];
    double ljk = originalEdgeLengths[he.next().edge()];
    double lki = originalEdgeLengths[he.next().next().edge()];
    double lil = originalEdgeLengths[he.twin().next().edge()];
    double llj = originalEdgeLengths[he.twin().next().next().edge()];
    return {lij, ljk, lki, lil, llj};
}

std::array<double, 5>
Triangulation::originalNeighborhoodEdgeLengths(Edge e) const {
    return originalNeighborhoodEdgeLengths(e.halfedge());
}

std::array<size_t, 5>
Triangulation::neighborhoodNormalCoordinates(Halfedge he) const {
    // Gather normal coordinates
    size_t nij = normalCoordinates[he.edge()];
    size_t njk = normalCoordinates[he.next().edge()];
    size_t nki = normalCoordinates[he.next().next().edge()];
    size_t nil = normalCoordinates[he.twin().next().edge()];
    size_t nlj = normalCoordinates[he.twin().next().next().edge()];
    return {nij, njk, nki, nil, nlj};
}

std::array<size_t, 5>
Triangulation::neighborhoodNormalCoordinates(Edge e) const {
    return neighborhoodNormalCoordinates(e.halfedge());
}
} // namespace CEPS
