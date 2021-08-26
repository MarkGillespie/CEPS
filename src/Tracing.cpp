#include "Tracing.h"

namespace CEPS {
MeshPath traceEdge(Edge e1, Triangulation& T1, Triangulation& T2) {
    auto vertexHalfedge = [&](Vertex v, size_t iH) {
        Halfedge he = v.halfedge();

        // Iterate counterclockwise
        for (size_t i = 0; i < iH; ++i) {
            he = he.next().next().twin();
        }
        return he;
    };

    auto directPath = [&](Halfedge he) {
        MeshPath path;
        path.points = {};
        path.start  = CornerPt{he.corner(), 1};
        path.end    = CornerPt{he.next().corner(), 1};
        return path;
    };

    // Use cached vertex indices
    T1.geo->requireVertexIndices();
    const VertexData<size_t>& vIdx = T1.geo->vertexIndices;
    T1.geo->unrequireVertexIndices();

    Vertex vTrace = src(e1);
    Vertex v      = T2.mesh->vertex(vIdx[vTrace]);

    Halfedge he = v.halfedge();
    do {
        size_t iStart = T2.roundaboutIndices[he];

        size_t em       = positivePart(T2.normalCoordinates[he.next().edge()] -
                                 T2.normalCoordinates[he.edge()] -
                                 T2.normalCoordinates[he.next().next().edge()]);
        size_t startInd = (T2.normalCoordinates[he.edge()] == 0) ? 1 : 0;
        size_t endInd =
            (T2.normalCoordinates[he.next().next().edge()] == 0) ? 1 : 0;

        size_t width = em + startInd + endInd;

        for (size_t iH = 0; iH < width; ++iH) {
            Halfedge heTrace = vertexHalfedge(vTrace, iStart + iH);
            if (heTrace != e1.halfedge()) continue;
            if (iH == width - 1 &&
                T2.normalCoordinates[he.next().next().edge()] == 0) {

                Halfedge pathHe = he.next().next().twin();
                return directPath(pathHe);

            } else if (iH == 0 && T2.normalCoordinates[he.edge()] == 0) {
                return directPath(he);
            } else {
                size_t idx = iH;

                if (T2.normalCoordinates[he.edge()] == 0) {
                    idx = iH - 1;
                }

                // traceTopologicalCurve takes 1-indexed indices
                return ImplementationDetails::traceTopologicalCurve(
                    he.next(), T2.normalCoordinates[he.edge()] + idx + 1,
                    T2.normalCoordinates);
            }
        }
        // orbit counterclockwise
        he = he.next().next().twin();
    } while (he != v.halfedge());


    std::cerr << "Something somewhere went horribly wrong" << vendl;
    he = v.halfedge();
    do {
        size_t iStart = T2.roundaboutIndices[he];

        size_t em       = positivePart(T2.normalCoordinates[he.next().edge()] -
                                 T2.normalCoordinates[he.edge()] -
                                 T2.normalCoordinates[he.next().next().edge()]);
        size_t startInd = (T2.normalCoordinates[he.edge()] == 0) ? 1 : 0;
        size_t endInd =
            (T2.normalCoordinates[he.next().next().edge()] == 0) ? 1 : 0;

        size_t width = em + startInd + endInd;
        std::cerr << "iStart: " << iStart << "\tem: " << em
                  << "\tstartInd: " << startInd << "\t endInd: " << endInd
                  << vendl;
        for (size_t iH = 0; iH < width; ++iH) {
            std::cerr << "\ttrying halfedge " << iStart + iH << " of "
                      << vTrace.degree() << " (aka " << T2.roundaboutDegrees[v]
                      << ")" << vendl;
        }
        // orbit counterclockwise
        he = he.next().next().twin();
    } while (he != v.halfedge());

    return MeshPath();
}

// Trace T1 over T2
EdgeData<MeshPath> traceTopologicalTriangulation(Triangulation& T1,
                                                 Triangulation& T2,
                                                 bool verbose) {
    EdgeData<MeshPath> tracedEdges(*T1.mesh);

    size_t N  = T2.mesh->nEdges();
    size_t iE = 0;
    if (verbose) std::cout << "\t\ttracing edges | 0 / " << N;


    for (Edge e : T1.mesh->edges()) {
        tracedEdges[e] = traceEdge(e, T1, T2);
        if (verbose)
            std::cout << "\r\t\ttracing edges | " << ++iE << " / " << N
                      << std::flush;
    }
    if (verbose) std::cout << std::endl;

    return tracedEdges;
}

// Trace T1 over T2
EdgeData<MeshPath> traceGeodesicTriangulation(Triangulation& T1,
                                              Triangulation& T2,
                                              GeometryType gType,
                                              bool verbose) {
    EdgeData<MeshPath> paths = traceTopologicalTriangulation(T1, T2, verbose);

    size_t N  = T2.mesh->nEdges();
    size_t iE = 0;
    if (verbose) std::cout << "\t\tstraightening edges | 0 / " << N;

    switch (gType) {
    case GeometryType::EUCLIDEAN:
        for (Edge e : T1.mesh->edges()) {
            straightenEuclidean(paths[e], T2);
            if (verbose)
                std::cout << "\r\t\tstraightening edges | " << ++iE << " / "
                          << N << std::flush;
        }
        break;
    case GeometryType::HYPERBOLIC:
        VertexData<double> negativeScaleFactors = -T2.logScaleFactors;
        for (Edge e : T1.mesh->edges()) {
            bool straighteningSucceeded =
                straightenHyperbolic(paths[e], T2, negativeScaleFactors);
            if (!straighteningSucceeded)
                straightenHyperbolicIteratively(paths[e], T2,
                                                negativeScaleFactors);
            if (verbose)
                std::cout << "\r\t\tstraightening edges | " << ++iE << " / "
                          << N << std::flush;
        }
        break;
    }
    if (verbose) std::cout << std::endl;

    return paths;
}

// Trace T1 over T2, doing as much computation on T1 as possible
// This turns out to be an easier problem in many instances
EdgeData<MeshPath> traceTransposedGeodesicTriangulation(Triangulation& T1,
                                                        Triangulation& T2,
                                                        GeometryType gType,
                                                        bool verbose) {
    using ImplementationDetails::transpose;

    EdgeData<MeshPath> pathsOnT2 =
        traceTopologicalTriangulation(T1, T2, verbose);

    EdgeData<MeshPath> pathsOnT1 = transpose(pathsOnT2, T1, T2);

    size_t N  = T2.mesh->nEdges();
    size_t iE = 0;
    if (verbose) std::cout << "\t\tstraightening edges | 0 / " << N;

    switch (gType) {
    case GeometryType::EUCLIDEAN:
        for (Edge e : T2.mesh->edges()) {
            straightenEuclidean(pathsOnT1[e], T1);
            if (verbose)
                std::cout << "\r\t\tstraightening edges | " << ++iE << " / "
                          << N << std::flush;
        }
        break;
    case GeometryType::HYPERBOLIC:
        VertexData<double> scaleFactors =
            T2.logScaleFactors.reinterpretTo(*T1.mesh);
        for (Edge e : T2.mesh->edges()) {
            bool straighteningSucceeded =
                straightenHyperbolic(pathsOnT1[e], T1, scaleFactors);
            if (!straighteningSucceeded)
                straightenHyperbolicIteratively(pathsOnT1[e], T1, scaleFactors);
            if (verbose)
                std::cout << "\r\t\tstraightening edges | " << ++iE << " / "
                          << N << std::flush;
        }
        break;
    }
    if (verbose) std::cout << std::endl;

    return transpose(pathsOnT1, T2, T1);
}

void straightenEuclidean(MeshPath& curve, Triangulation& T) {
    if (curve.points.empty()) {
        return;
    }

    T.geo->requireEdgeLengths();
    T.geo->requireCornerAngles();

    auto layOutFirstFace = [&](Halfedge he) {
        double len   = T.geo->edgeLengths[he.edge()];
        Vector2 pos1 = Vector2{0, 0};
        Vector2 pos2 = Vector2{len, 0};
        len          = T.geo->edgeLengths[he.next().next().edge()];
        double angle = T.geo->cornerAngle(he.corner());
        Vector2 pos0 = Vector2{len * cos(angle), len * sin(angle)};

        return std::make_tuple(pos0, pos1, pos2);
    };

    auto nextAngle = [&](Halfedge he, Halfedge prevHe, double prevAngle) {
        double angle = T.geo->cornerAngle(prevHe.twin().next().corner());
        double theta = prevAngle - angle;

        if (prevHe.twin().next() == he) {
            return std::make_tuple(theta, theta);
        } else if (prevHe.twin().next().next() == he) {
            double nextAngle =
                T.geo->cornerAngle(prevHe.twin().next().next().corner());
            return std::make_tuple(theta, theta + M_PI - nextAngle);
        } else {
            throw_verbose_runtime_error(
                "Path connectivity error in "
                "Tracing.cpp:straightenEuclidean:nextAngle");
        }
    };

    auto layOutNextFaceFromAngle =
        [&](Halfedge he, Halfedge prevHe,
            const std::array<Vector2, 3>& prevPositions, double theta) {
            double len = T.geo->edgeLengths[prevHe.twin().next().edge()];

            Vector2 newPos =
                prevPositions[1] + Vector2{len * cos(theta), len * sin(theta)};

            if (prevHe.twin().next() == he) {
                return std::tuple<Vector2, Vector2, Vector2>{
                    prevPositions[2], prevPositions[1], newPos};
            } else if (prevHe.twin().next().next() == he) {
                return std::tuple<Vector2, Vector2, Vector2>{
                    prevPositions[1], newPos, prevPositions[2]};
            } else {
                throw_verbose_runtime_error(
                    "Path connectivity error in "
                    "Tracing.cpp:straightenEuclidean:layOutNextFaceFromAngle");
            }
        };

    // TODO: replace with placeTrianglePoint or whatever it's called
    auto layOutNextFace = [&](Halfedge he, Halfedge prevHe,
                              const std::array<Vector2, 3>& prevPositions) {
        double angle = T.geo->cornerAngle(prevHe.twin().next().corner());
        double len   = T.geo->edgeLengths[prevHe.twin().next().edge()];
        Vector2 disp = prevPositions[2] - prevPositions[1];
        double theta = atan2(disp.y, disp.x) - angle;

        Vector2 newPos =
            prevPositions[1] + Vector2{len * cos(theta), len * sin(theta)};


        // TODO: it's kind of wasteful to return a whole tuple, but it makes
        // my life easier for now
        if (prevHe.twin().next() == he) {
            return std::tuple<Vector2, Vector2, Vector2>{
                prevPositions[2], prevPositions[1], newPos};
        } else if (prevHe.twin().next().next() == he) {
            return std::tuple<Vector2, Vector2, Vector2>{
                prevPositions[1], newPos, prevPositions[2]};
        } else {
            exit(1);
        }
    };

    std::vector<std::array<Vector2, 2>> edgePositions;
    std::vector<double> angles;

    Vector2 base, src, dst;

    Halfedge prevHe = curve.points[0].halfedge;

    std::tie(base, src, dst) = layOutFirstFace(prevHe);
    edgePositions.push_back({src, dst});

    Vector2 start = base;

    double heAngle = 0, refAngle = 0;
    for (size_t iC = 1; iC < curve.points.size(); ++iC) {
        Halfedge he                 = curve.points[iC].halfedge;
        std::tie(refAngle, heAngle) = nextAngle(he, prevHe, heAngle);

        angles.push_back(refAngle);
        prevHe = he;
    }

    prevHe = curve.points[0].halfedge;
    for (size_t iC = 1; iC < curve.points.size(); ++iC) {
        Halfedge he              = curve.points[iC].halfedge;
        std::tie(base, src, dst) = layOutNextFaceFromAngle(
            he, prevHe, {base, src, dst}, angles[iC - 1]);

        edgePositions.push_back({src, dst});
        prevHe = he;
    }

    Halfedge almostLastHe = curve.points[curve.points.size() - 1].halfedge;
    // lastHe points to the ending vertex
    Halfedge lastHe = almostLastHe.twin().next();
    std::tie(base, src, dst) =
        layOutNextFace(lastHe, almostLastHe, {base, src, dst});
    Vector2 end = dst;

    bool displayed  = false;
    double curveLen = (start - end).norm();
    curve.baryCoords.clear();
    curve.baryCoords.reserve(curve.points.size());
    for (size_t iC = 0; iC < curve.points.size(); ++iC) {
        Halfedge he = curve.points[iC].halfedge;

        double t = intersectionTime(edgePositions[iC][0], edgePositions[iC][1],
                                    start, end);

        // Clamp t to [0,1]
        t = fmin(fmax(t, 0), 1);

        curve.points[iC] = HalfedgePt{he, bary2(t)};

        Vector2 intersection =
            t * edgePositions[iC][0] + (1 - t) * edgePositions[iC][1];

        double r = (intersection - start).norm();
        curve.baryCoords.push_back(r / curveLen);
    }
}

// Straighten Hyperbolic geodesic over T, according to given scale factors
// return false if straightening fails
bool straightenHyperbolic(MeshPath& curve, Triangulation& T,
                          const VertexData<double>& logScaleFactors) {

    if (curve.points.empty()) {
        Vertex vStart = curve.start.corner.vertex();
        Vertex vEnd   = curve.end.corner.vertex();

        curve.start.bary = exp(logScaleFactors[vStart]);
        curve.end.bary   = exp(logScaleFactors[vEnd]);

        return true;
    } else if (curve.points.size() == 1) {

        Vertex vStart = curve.start.corner.vertex();
        Vertex vEnd   = curve.end.corner.vertex();

        Halfedge he = curve.points[0].halfedge;

        auto len = [&](Halfedge he) -> double {
            return T.geo->inputEdgeLengths[he.edge()];
        };

        double startScale = exp(logScaleFactors[vStart]);
        double endScale   = exp(logScaleFactors[vEnd]);

        double l_ab    = len(he);
        double l_bs    = len(he.next()) * sqrt(startScale);
        double l_sa    = len(he.next().next()) * sqrt(startScale);
        double l_ae    = len(he.twin().next()) * sqrt(endScale);
        double l_eb    = len(he.twin().next().next()) * sqrt(endScale);
        double l_se    = (l_bs * l_ae + l_eb * l_sa) / l_ab;
        double det_sab = l_ab * l_bs * l_sa / 2.;
        double det_eba = l_eb * l_ab * l_ae / 2.;
        double det_bse = l_bs * l_se * l_eb / 2.;
        double det_aes = l_ae * l_se * l_sa / 2.;
        double dot_se  = -0.5 * pow(l_se, 2);
        double dot_sa  = -0.5 * pow(l_sa, 2);
        double dot_bs  = -0.5 * pow(l_bs, 2);

        double pathBaryCoord = det_sab / (det_sab + det_eba);

        Vector2 hedgeBary{det_bse, det_aes};
        hedgeBary    = normalizeBary(hedgeBary);
        double scale = (pathBaryCoord * dot_se) /
                       (hedgeBary.x * dot_sa + hedgeBary.y * dot_bs);
        hedgeBary *= scale;

        curve.baryCoords = {pathBaryCoord};

        if (l_se != l_se || pathBaryCoord != pathBaryCoord ||
            hedgeBary.x != hedgeBary.x || hedgeBary.y != hedgeBary.y)
            return false;

        curve.start.bary = startScale;
        curve.end.bary   = endScale;
        curve.points[0]  = HalfedgePt{he, hedgeBary};

        return true;
    }

    auto layOutHyperboloidTri =
        [&](Halfedge he) -> std::tuple<Vector3, Vector3, Vector3> {
        double u = T.geo->inputEdgeLengths[he.next().next().edge()];
        double v = T.geo->inputEdgeLengths[he.next().edge()];
        double d = T.geo->inputEdgeLengths[he.edge()];

        return ImplementationDetails::computeLightConeCoords(v, d, u);
    };

    auto layOutNextHyperboloidTri =
        [&](Halfedge he, Halfedge prevHe,
            const std::array<Vector3, 3>& prevVertices)
        -> std::tuple<Vector3, Vector3, Vector3> {
        double u = T.geo->inputEdgeLengths[prevHe.next().next().edge()];
        double v = T.geo->inputEdgeLengths[prevHe.next().edge()];
        double d = T.geo->inputEdgeLengths[prevHe.edge()];
        double x = T.geo->inputEdgeLengths[prevHe.twin().next().edge()];
        double y = T.geo->inputEdgeLengths[prevHe.twin().next().next().edge()];
        double diag = (u * y + v * x) / d;

        Vector3 newPt = ImplementationDetails::placeFourthLightConePoint(
            prevVertices[2], prevVertices[0], prevVertices[1], y, diag, x);

        auto prevDistErr = [&](double dist, size_t vi, size_t vj) {
            return abs(dist - ImplementationDetails::lorentzDist(
                                  prevVertices[vi], prevVertices[vj]));
        };

        if (prevDistErr(u, 0, 1) > 1e-4 || prevDistErr(v, 0, 2) > 1e-4 ||
            prevDistErr(d, 1, 2) > 1e-4) {
            Vector3 inf{std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};

            return {inf, inf, inf};
        }

        if (prevHe.twin().next() == he) {
            return {prevVertices[2], prevVertices[1], newPt};
        } else if (prevHe.twin().next().next() == he) {
            return {prevVertices[1], newPt, prevVertices[2]};
        } else {
            Vector3 inf{std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
            std::cerr << "bad connectivity!" << vendl;
            return {inf, inf, inf};
        }
    };

    std::vector<std::array<Vector3, 2>> edgePositions;

    Vector3 base, src, dst;

    Halfedge prevHe = curve.points[0].halfedge;

    std::tie(base, src, dst) = layOutHyperboloidTri(prevHe);
    edgePositions.push_back({src, dst});

    double startScale   = exp(logScaleFactors[curve.start.corner.vertex()]);
    Vector3 start       = base;
    Vector3 scaledStart = startScale * start;

    for (size_t iP = 0; iP + 1 < curve.points.size(); ++iP) {
        Halfedge he     = curve.points[iP].halfedge;
        Halfedge heNext = curve.points[iP + 1].halfedge;
        if (heNext != he.twin().next() && heNext != he.twin().next().next()) {
            std::cerr << "Tracing connectivity error?" << vendl;
        }
    }

    for (size_t iC = 1; iC < curve.points.size(); ++iC) {
        Halfedge he = curve.points[iC].halfedge;
        std::tie(base, src, dst) =
            layOutNextHyperboloidTri(he, prevHe, {base, src, dst});
        if (!std::isfinite(base.norm()) || !std::isfinite(src.norm()) ||
            !std::isfinite(dst.norm())) {
            // std::cerr << "\t err in middle! [" << iC << "]" << vendl;
            return false;
        }
        edgePositions.push_back({src, dst});
        prevHe = he;
    }

    Halfedge almostLastHe = curve.points[curve.points.size() - 1].halfedge;
    // lastHe points to the ending vertex
    Halfedge lastHe = almostLastHe.twin().next();
    std::tie(base, src, dst) =
        layOutNextHyperboloidTri(lastHe, almostLastHe, {base, src, dst});
    if (!std::isfinite(base.norm()) || !std::isfinite(src.norm()) ||
        !std::isfinite(dst.norm())) {
        // std::cerr << "\t err at end!" << vendl;
        return false;
    }
    edgePositions.push_back({src, dst});
    Vector3 end       = dst;
    double endScale   = exp(logScaleFactors[curve.end.corner.vertex()]);
    Vector3 scaledEnd = endScale * end;

    curve.baryCoords.clear();
    curve.baryCoords.reserve(curve.points.size());

    curve.start.bary = startScale;
    curve.end.bary   = endScale;

    for (size_t iC = 0; iC < curve.points.size(); ++iC) {
        Halfedge he = curve.points[iC].halfedge;

        double uSrc = logScaleFactors[he.vertex()];
        double uDst = logScaleFactors[he.twin().vertex()];

        Vector2 hIntersectionBary =
            ImplementationDetails::homogeneousProjectiveIntersectionTime(
                edgePositions[iC][0], edgePositions[iC][1], scaledStart,
                scaledEnd);
        Vector2 intersectionBary = normalizeBary(hIntersectionBary);
        double t                 = intersectionBary[0];
        double s = ImplementationDetails::projectiveIntersectionTime(
            scaledStart, scaledEnd, edgePositions[iC][0], edgePositions[iC][1]);

        if (t != t || t < 0 || t > 1) {
            return false;
        }

        curve.points[iC] = HalfedgePt{he, hIntersectionBary};

        Vector3 intersection =
            t * edgePositions[iC][0] + (1 - t) * edgePositions[iC][1];
        curve.baryCoords.push_back(1 - s);

        if (s < 0 || s > 1) return false;
    }
    return true;
}

bool straightenHyperbolicIteratively(MeshPath& curve, Triangulation& T,
                                     const VertexData<double>& logScaleFactors,
                                     double tol) {
    using ImplementationDetails::twin;

    T.geo->requireEdgeLengths();

    Vertex vS = curve.start.corner.vertex();
    Vertex vE = curve.end.corner.vertex();

    double srcScale = exp(logScaleFactors[vS]);
    double dstScale = exp(logScaleFactors[vE]);

    curve.start.bary = srcScale;
    curve.end.bary   = dstScale;

    // Exit early if there are no interior points to straighten
    if (curve.points.empty()) return true;

    double maxChange = 2 * tol;
    size_t iter      = 0;

    curve.baryCoords.clear();
    for (size_t iP = 0; iP < curve.points.size(); ++iP) {
        curve.baryCoords.push_back((double)(iP + 1) /
                                   (double)(curve.points.size() + 1));
    }

    // Components in order base, a, b, op
    Eigen::Vector4d startPoint, endPoint;
    Eigen::Vector4d base(1, 0, 0, 0);
    Eigen::Vector4d a(0, 1, 0, 0);
    Eigen::Vector4d b(0, 0, 1, 0);
    Eigen::Vector4d opp(0, 0, 0, 1);

    const EdgeData<double>& len = T.geo->edgeLengths;

    // Temporary variables
    std::vector<std::vector<double>> lengths(4, std::vector<double>(4, 0));
    std::array<std::array<std::array<double, 4>, 4>, 4> basisAbsDet;


    while (maxChange > tol &&
           iter < 1e10 * curve.points.size() * curve.points.size()) {
        maxChange = -1;
        iter++;
        for (size_t i = 0; i < curve.points.size(); ++i) {
            HalfedgePt prevPt, pt, nextPt;
            pt = curve.points[i];
            if (i > 0) {
                prevPt = twin(curve.points[i - 1]);
            }

            if (i + 1 < curve.points.size()) {
                nextPt = curve.points[i + 1];
            }

            double l_base_a  = len[pt.halfedge.next().next().edge()];
            double l_base_b  = len[pt.halfedge.next().edge()];
            double l_a_b     = len[pt.halfedge.edge()];
            double l_op_a    = len[pt.halfedge.twin().next().edge()];
            double l_op_b    = len[pt.halfedge.twin().next().next().edge()];
            double l_base_op = (l_base_a * l_op_b + l_base_b * l_op_a) / l_a_b;

            lengths[0][1] = l_base_a;
            lengths[0][2] = l_base_b;
            lengths[0][3] = l_base_op;

            lengths[1][0] = l_base_a;
            lengths[1][2] = l_a_b;
            lengths[1][3] = l_op_a;

            lengths[2][0] = l_base_b;
            lengths[2][1] = l_a_b;
            lengths[2][3] = l_op_b;

            lengths[3][0] = l_base_op;
            lengths[3][1] = l_op_a;
            lengths[3][2] = l_op_b;

            for (size_t iX = 0; iX < 4; ++iX) {
                for (size_t iY = 0; iY < 4; ++iY) {
                    for (size_t iZ = 0; iZ < 4; ++iZ) {
                        basisAbsDet[iX][iY][iZ] =
                            abs(0.5 * lengths[iX][iY] * lengths[iY][iZ] *
                                lengths[iZ][iX]);
                    }
                }
            }

            auto absDet = [&](const Eigen::Vector4d& x,
                              const Eigen::Vector4d& y,
                              const Eigen::Vector4d& z) -> double {
                double total = 0;
                for (size_t iX = 0; iX < 4; ++iX) {
                    for (size_t iY = 0; iY < 4; ++iY) {
                        for (size_t iZ = 0; iZ < 4; ++iZ) {
                            double coeff = x(iX) * y(iY) * z(iZ);
                            total += coeff * basisAbsDet[iX][iY][iZ];
                        }
                    }
                }
                return total;
            };

            auto basisDot = [&](size_t iX, size_t iY) -> double {
                return -0.5 * pow(lengths[iX][iY], 2);
            };
            auto dot = [&](const Eigen::Vector4d& x,
                           const Eigen::Vector4d& y) -> double {
                double total = 0;
                for (size_t iX = 0; iX < 4; ++iX) {
                    for (size_t iY = 0; iY < 4; ++iY) {
                        double coeff = x(iX) * y(iY);
                        total += coeff * basisDot(iX, iY);
                    }
                }
                return total;
            };

            double startScale, endScale, startCurveBary, endCurveBary;

            if (i == 0) {
                startScale     = srcScale;
                startPoint     = Eigen::Vector4d(1, 0, 0, 0);
                startCurveBary = 0;
            } else if (prevPt.halfedge.next().vertex() ==
                       pt.halfedge.vertex()) {
                // start on edge u
                startScale     = prevPt.bary.x + prevPt.bary.y;
                startPoint     = Eigen::Vector4d(prevPt.bary.x / startScale,
                                             prevPt.bary.y / startScale, 0, 0);
                startCurveBary = curve.baryCoords[i - 1];
            } else {
                // start on edge v
                startScale     = prevPt.bary.x + prevPt.bary.y;
                startPoint     = Eigen::Vector4d(prevPt.bary.y / startScale, 0,
                                             prevPt.bary.x / startScale, 0);
                startCurveBary = curve.baryCoords[i - 1];
            }

            if (i + 1 == curve.points.size()) {
                endScale     = dstScale;
                endPoint     = Eigen::Vector4d(0, 0, 0, 1);
                endCurveBary = 1;
            } else if (nextPt.halfedge.vertex() == pt.halfedge.vertex()) {
                // end on edge x
                endScale     = nextPt.bary.x + nextPt.bary.y;
                endPoint     = Eigen::Vector4d(0, nextPt.bary.x / endScale, 0,
                                           nextPt.bary.y / endScale);
                endCurveBary = curve.baryCoords[i + 1];
            } else {
                // end on edge y
                endScale     = nextPt.bary.x + nextPt.bary.y;
                endPoint     = Eigen::Vector4d(0, 0, nextPt.bary.y / endScale,
                                           nextPt.bary.x / endScale);
                endCurveBary = curve.baryCoords[i + 1];
            }

            // Special case for a curve which just barely clips off a corner
            if (i > 0 && i + 1 < curve.points.size()) {
                if (prevPt.halfedge.next().vertex() ==
                    nextPt.halfedge.vertex()) {
                    if (prevPt.bary.y / prevPt.bary.x > 1e12 &&
                        nextPt.bary.x / nextPt.bary.y > 1e12) {

                        double avgScale = (prevPt.bary.x + prevPt.bary.y +
                                           nextPt.bary.x + nextPt.bary.y) /
                                          2.;

                        verbose_assert(!std::isnan(avgScale),
                                       "avgScale is nan?");

                        Vector2 newBary{avgScale, 0};

                        curve.points[i] = HalfedgePt{pt.halfedge, newBary};
                        curve.baryCoords[i] =
                            (startCurveBary + endCurveBary) / 2.;

                        continue;
                    }
                } else if (nextPt.halfedge.next().vertex() ==
                           prevPt.halfedge.vertex()) {

                    if (prevPt.bary.y / prevPt.bary.x < 1e-12 &&
                        nextPt.bary.x / nextPt.bary.y < 1e-12) {
                        double avgScale = (prevPt.bary.x + prevPt.bary.y +
                                           nextPt.bary.x + nextPt.bary.y) /
                                          2.;
                        verbose_assert(!std::isnan(avgScale),
                                       "avgScale is nan?");

                        Vector2 newBary{0, avgScale};

                        curve.points[i] = HalfedgePt{pt.halfedge, newBary};
                        curve.baryCoords[i] =
                            (startCurveBary + endCurveBary) / 2.;

                        continue;
                    }
                }
            }

            Vector2 oldBary = curve.points[i].bary;

            Vector2 intrinsicBary{absDet(b, startPoint, endPoint),
                                  absDet(a, startPoint, endPoint)};
            // For numerical stability, we normalize intrinsicBary
            intrinsicBary = normalizeBary(intrinsicBary);

            // Solve for scale factor at intersection, and the coefficients
            // describing the intersection as a convex combination of
            // startPoint and endPoint
            double g1 = endScale * dot(endPoint, b);
            double g2 = endScale * dot(endPoint, a);

            double a1 = intrinsicBary.x * dot(a, b);
            double a2 = intrinsicBary.y * dot(a, b);

            double b1 =
                endScale * dot(endPoint, b) - startScale * dot(startPoint, b);
            double b2 =
                endScale * dot(endPoint, a) - startScale * dot(startPoint, a);

            // If system is degenerate, obtain second equation by dotting with
            // startPoint instead of dotting with a
            if (abs(a1 * b2 - a2 * b1) < 1e-4) {
                a2 = intrinsicBary.x * startScale * dot(a, startPoint) +
                     intrinsicBary.y * startScale * dot(b, startPoint);
                b2 = endScale * startScale * dot(endPoint, startPoint) -
                     startScale * startScale * dot(startPoint, startPoint);
                g2 = endScale * startScale * dot(endPoint, startPoint);
            }

            double intrinsicScale = (g1 * b2 - g2 * b1) / (a1 * b2 - a2 * b1);
            double sigma          = (a1 * g2 - a2 * g1) / (a1 * b2 - a2 * b1);


            Vector2 newBary = abs(intrinsicScale) * intrinsicBary;

            double newPathBary =
                sigma * startCurveBary + (1 - sigma) * endCurveBary;

            if (newBary.x < 0 || newBary.y < 0) {
                std::cout << "invalid barycentric coordinate " << newBary
                          << vendl;
            }
            newBary.x = fmax(0, newBary.x);
            newBary.y = fmax(0, newBary.y);

            double baryChange =
                abs(oldBary.x - newBary.x) + abs(oldBary.y - newBary.y);
            double pathBaryChange = abs(newPathBary - curve.baryCoords[i]);
            maxChange = fmax(maxChange, fmax(baryChange, pathBaryChange));

            curve.points[i]     = HalfedgePt{pt.halfedge, newBary};
            curve.baryCoords[i] = newPathBary;
        }
    }
    // std::cerr << "INTRINSIC STRAIGHTENING FINISHED AFTER " << iter
    //           << " ITERATIONS" << vendl;
    return maxChange > tol;
}

namespace ImplementationDetails {

// Flip the MeshPath to go in the other direction
MeshPath twin(MeshPath path) {

    // Reverse start and end
    CornerPt oldStart = path.start;
    CornerPt oldEnd   = path.end;
    path.start        = oldEnd;
    path.end          = oldStart;

    // Reverse point order
    std::reverse(path.points.begin(), path.points.end());

    if (!path.points.empty()) {
        // Flip HalfedgePts
        for (size_t i = 0; i < path.points.size(); ++i) {
            path.points[i] = twin(path.points[i]);
        }
    }

    if (!path.baryCoords.empty()) {
        size_t N = path.baryCoords.size();
        std::vector<double> new_baryCoords;
        for (size_t i = 0; i < N; ++i) {
            new_baryCoords.push_back(1 - path.baryCoords[N - i - 1]);
        }
        path.baryCoords = new_baryCoords;
    }

    return path;
}

// Return an equivalent HalfedgePt where he = he.edge().halfedge()
HalfedgePt canonicalize(HalfedgePt hPt) {
    if (hPt.halfedge.edge().halfedge() == hPt.halfedge) {
        return hPt;
    } else {
        return twin(hPt);
    }
}
HalfedgePt twin(HalfedgePt hPt) {
    return HalfedgePt{hPt.halfedge.twin(), Vector2{hPt.bary.y, hPt.bary.x}};
}

EdgeData<MeshPath> transpose(const EdgeData<MeshPath>& paths, Triangulation& T1,
                             Triangulation& T2) {

    EdgeData<MeshPath> t2overT1(*T2.mesh);

    EdgeData<std::vector<std::tuple<HalfedgePt, Halfedge, double>>>
        t2EdgePoints(*T2.mesh);

    VertexData<size_t> vIdx2 = T2.mesh->getVertexIndices();
    VertexData<size_t> vIdx1 = T1.mesh->getVertexIndices();
    std::vector<double> expU(T2.mesh->nVertices(), 1);
    std::unordered_map<Edge, Halfedge> t2toT1;

    for (Edge e : T1.mesh->edges()) {
        const MeshPath& path = paths[e];
        if (path.points.empty()) {
            Halfedge he2 = path.start.corner.halfedge();
            if (he2 == Halfedge()) {
                if (path.start.corner == Corner())
                    std::cout << "uninitialized corner" << vendl;
                throw_verbose_runtime_error(
                    "tried to transpose uninitialized halfedge");
            }
            if (he2 == he2.edge().halfedge()) {
                t2toT1[he2.edge()] = e.halfedge();
            } else {
                t2toT1[he2.edge()] = e.halfedge().twin();
            }
        } else {
            for (size_t iP = 0; iP < path.points.size(); ++iP) {
                HalfedgePt hePt     = path.points[iP];
                Halfedge crossingHe = hePt.halfedge;
                if (crossingHe == crossingHe.edge().halfedge()) {
                    t2EdgePoints[crossingHe.edge()].push_back(
                        std::make_tuple(canonicalize(hePt), e.halfedge().twin(),
                                        1 - path.baryCoords[iP]));
                } else {
                    t2EdgePoints[crossingHe.edge()].push_back(std::make_tuple(
                        canonicalize(hePt), e.halfedge(), path.baryCoords[iP]));
                }
            }
        }
        expU[vIdx2[path.start.corner.vertex()]] = 1. / path.start.bary;
        expU[vIdx2[path.end.corner.vertex()]]   = 1. / path.end.bary;
    }

    for (Edge e : T2.mesh->edges()) {
        if (t2EdgePoints[e].empty()) {
            Halfedge he1 = t2toT1[e];

            Corner start = he1.corner();
            Corner end   = he1.next().corner();

            t2overT1[e].start = CornerPt{start, expU[vIdx1[start.vertex()]]};
            t2overT1[e].end   = CornerPt{end, expU[vIdx1[end.vertex()]]};

            t2overT1[e].baryCoords = {};
        } else {
            auto cmpBary =
                [](const std::tuple<HalfedgePt, Halfedge, double>& a,
                   const std::tuple<HalfedgePt, Halfedge, double>& b) -> bool {
                Vector2 aBary = normalizeBary(std::get<0>(a).bary);
                Vector2 bBary = normalizeBary(std::get<0>(b).bary);

                return aBary.y < bBary.y;
            };

            std::sort(std::begin(t2EdgePoints[e]), std::end(t2EdgePoints[e]),
                      cmpBary);

            Corner start =
                std::get<1>(t2EdgePoints[e][0]).next().next().corner();
            size_t n   = t2EdgePoints[e].size();
            Corner end = std::get<1>(t2EdgePoints[e][n - 1])
                             .twin()
                             .next()
                             .next()
                             .corner();
            t2overT1[e].start = CornerPt{start, expU[vIdx1[start.vertex()]]};
            t2overT1[e].end   = CornerPt{end, expU[vIdx1[end.vertex()]]};

            HalfedgePt hPt;
            Halfedge overlayHe;
            double overlayBary;
            for (const auto& p : t2EdgePoints[e]) {
                std::tie(hPt, overlayHe, overlayBary) = p;
                double totalScale                     = hPt.bary.x + hPt.bary.y;

                t2overT1[e].points.push_back(HalfedgePt{
                    overlayHe,
                    Vector2{1 - overlayBary, overlayBary} / totalScale});
                t2overT1[e].baryCoords.push_back(hPt.bary.y / totalScale);
            }
        }
    }

    return t2overT1;
}

bool normalCoordinateTriangleInequalityViolation(Face f, Halfedge& violatingHe,
                                                 const EdgeData<size_t>& n) {
    Halfedge ha = f.halfedge();
    Halfedge hb = ha.next();
    Halfedge hc = hb.next();

    size_t a = n[ha.edge()];
    size_t b = n[hb.edge()];
    size_t c = n[hc.edge()];

    if (a > b + c) {
        violatingHe = ha;
        return true;
    } else if (b > c + a) {
        violatingHe = hb;
        return true;
    } else if (c > a + b) {
        violatingHe = hc;
        return true;
    } else {
        return false;
    }
}

// He is the twin of the halfedge in a face that the curve starts from,
// index is the index of the curve, with 1 being closest to the source
// vertex of he. The returned path starts on He, and ends at some vertex
// WARNING: 1-indexed
MeshPath traceTopologicalCurve(Halfedge he, size_t index,
                               const EdgeData<size_t>& n) {
    MeshPath path;

    path.start = CornerPt{he.next().next().corner(), 1};
    do {
        double t = (double)index / (double)(n[he.edge()] + 1);
        path.points.push_back(HalfedgePt{he, Vector2{1 - t, t}});
        std::tie(he, index) = nextEdge(he, index, n);
    } while (index > 0);

    path.end = CornerPt{he.next().next().corner(), 1};

    double N = path.points.size() + 1;
    path.baryCoords.reserve(path.points.size());
    for (size_t iP = 0; iP < path.points.size(); ++iP)
        path.baryCoords.push_back((double)(iP + 1) / N);

    return path;
}

std::tuple<Halfedge, size_t> nextEdge(Halfedge he, size_t p,
                                      const EdgeData<size_t>& n) {
    he         = he.twin();
    size_t nkj = n[he.edge()];
    size_t njl = n[he.next().edge()];
    size_t nlk = n[he.next().next().edge()];

    Halfedge violatingHe;
    if (!normalCoordinateTriangleInequalityViolation(he.face(), violatingHe,
                                                     n)) {
        // Case I
        size_t ml = (njl + nlk - nkj) / 2;
        if (p <= njl - ml) {
            return std::make_tuple(he.next(), p);
        } else {
            size_t mj = (njl + nkj - nlk) / 2;
            return std::make_tuple(he.next().next(), p - mj + ml);
        }
    } else if (violatingHe == he) {
        // Case II
        if (p <= njl) {
            return std::make_tuple(he.next(), p);
        } else if (p <= nkj - nlk) {
            return std::make_tuple(he, 0);
        } else {
            return std::make_tuple(he.next().next(), p - (nkj - nlk));
        }
    } else {
        // Case III
        if (violatingHe == he.next()) {
            return std::make_tuple(violatingHe, p);
        } else {
            return std::make_tuple(violatingHe, p + nlk - nkj);
        }
    }
    return std::make_tuple(he, p);
}

/* Inscribes a triangle with side lengths u, v, d into the
 * hyperboloid model. Returns vertex positions w0, w1, w2 on the light cone so
 * that |w1-w0| = d, |w2 - w1| = v and |w0 - w2| = u
 * (|.| denotes the Lorentz metric on R3, where the time coordinate is z)
 */
std::tuple<Vector3, Vector3, Vector3> computeLightConeCoords(double u, double v,
                                                             double d) {

    double d2 = pow(d, 2);
    double v2 = pow(v, 2);
    double u2 = pow(u, 2);

    double s1 = safeSqrt(d2 * u2 / (3 * v2));
    double s2 = safeSqrt(d2 * v2 / (3 * u2));
    double s3 = safeSqrt(u2 * v2 / (3 * d2));

    // Each disinct pair w_i, w_j has inner product -3/2. By rescaling each w_i
    // by s_i, we make these inner products (-3/2) * s_i * s_j instead. We set
    // s_i so that -3/2 s_i s_j = -1/2 length^2 for the length of the
    // appropriate edge
    Vector3 w1{1, 0, 1};
    Vector3 w2{cos(2 * M_PI / 3), sin(2 * M_PI / 3), 1};
    Vector3 w3{cos(2 * M_PI * 2 / 3), sin(2 * M_PI * 2 / 3), 1};


    return std::make_tuple(s1 * w1, s2 * w2, s3 * w3);
}


/* Computes the Lorentz metric on R3, IE
 * v.x * w.x + v.y * w.y - v.z * w.z
 */
double lorentz(Vector3 v, Vector3 w) {
    return v.x * w.x + v.y * w.y - v.z * w.z;
}

/* Returns lorentz(v, v)
 */
double lorentzNorm2(Vector3 v) { return lorentz(v, v); }

/* Computes the squared distance between v and w in the Lorentz metric
 */
double lorentzDist2(Vector3 v, Vector3 w) { return lorentzNorm2(v - w); }

/* Computes the distance between v and w in the Lorentz metric
 */
double lorentzDist(Vector3 v, Vector3 w) { return sqrt(lorentzDist2(v, w)); }

/* Takes in lightlike vectors v1, v2, v3 and distance l1, l2, l3. Computes a
 * point v4 on the light cone so that the (Minsk) distance to v1 is l1, etc
 */
Vector3 placeFourthLightConePoint(Vector3 v1, Vector3 v2, Vector3 v3, double l1,
                                  double l2, double l3) {
    // Compute pairwise distances
    double l14 = l1;
    double l24 = l2;
    double l34 = l3;
    double l12 = lorentzDist(v1, v2);
    double l13 = lorentzDist(v1, v3);
    double l23 = lorentzDist(v2, v3);

    // Coefficients such that \sum_i a_i v_i = 0
    double a1 = -l23 * l34 / l13;
    double a2 = l14 * l34 / l24;
    double a3 = -l14 * l12 / l13;
    double a4 = l23 * l12 / l24;

    Vector3 v4 = -(a1 * v1 + a2 * v2 + a3 * v3) / a4;

    return v4;
}

// Returns the intersection (v.x, v.y) of the projective line a1-a2 with
// projective line b1-b2 in barycentric coordinates along line a1-a2 (i.e.
// v.x * a1 + v.y * a2 lies on line b1-b2 as well). In this version, a1, a2, b1
// and b2 must be lightlike vectors
Vector2 homogeneousProjectiveIntersectionTime(Vector3 a1, Vector3 a2,
                                              Vector3 b1, Vector3 b2) {
    double sa = projectiveIntersectionTime(a1, a2, b1, b2);
    double sb = projectiveIntersectionTime(b1, b2, a1, a2);

    Vector3 intersectionA = sa * a1 + (1 - sa) * a2;
    Vector3 intersectionB = sb * b1 + (1 - sb) * b2;

    double scale = intersectionB.norm() / intersectionA.norm();

    return {scale * sa, scale * (1 - sa)};
}

// Returns the intersection of the projective line a1-a2 with projective line
// b1-b2 in barycentric coordinates along line a1-a2 (i.2. t*a1 + (1-t)a2 lies
// on line b1-b2 as well (up to scalar multiplication))
// a1, a2 must be lightlike vectors
double projectiveIntersectionTime(Vector3 a1, Vector3 a2, Vector3 b1,
                                  Vector3 b2) {
    auto det = [](Vector3 a, Vector3 b, Vector3 c) {
        return a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y -
               a.z * b.y * c.x - a.y * b.x * c.z - a.x * b.z * c.y;
    };

    double c1 = abs(det(a2, b2, b1));
    double c2 = abs(det(a1, b1, b2));

    return c1 / (c1 + c2);
}
} // namespace ImplementationDetails
} // namespace CEPS
