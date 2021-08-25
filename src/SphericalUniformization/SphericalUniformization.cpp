#include "SphericalUniformization.h"


namespace CEPS {


SphericalUniformizationResult
hemisphericalUniformize(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
                        bool viz) {
    // TODO: double mesh
    throw_verbose_runtime_error("Not implemented yet");
}

SphericalUniformizationResult sphericalUniformize(ManifoldSurfaceMesh& mesh,
                                                  VertexPositionGeometry& geo,
                                                  Vertex vInfinity, bool viz) {
    using ImplementationDetails::computeCommonRefinementAndMatrix;
    using ImplementationDetails::layOutTriangulation;
    using ImplementationDetails::mobiusCenter;

    // ==========================================
    //      Initialize all of the meshes
    // ==========================================
    PartiallyDecoratedTriangulation Ta(mesh, geo);
    VertexData<Vector3> vertexPositions =
        geo.inputVertexPositions.reinterpretTo(*Ta.mesh);

    // TODO: filter front faces
    std::set<Face> frontFaces;
    for (Face f : Ta.mesh->faces()) frontFaces.insert(f);

    // Copy Ta into Tb and then flip to intrinsic Delaunay
    PartiallyDecoratedTriangulation Tb(Ta, false);
    Tb.flipToDelaunay(GeometryType::EUCLIDEAN);

    // Copy Tb into Tc
    PartiallyDecoratedTriangulation Tc(Tb, true);

    if (vInfinity == Vertex()) {
        VectorHeatMethodSolver solver(*Tc.geo, 1);
        Tc.geo->requireVertexDualAreas();
        vInfinity =
            findCenter(*Tc.mesh, *Tc.geo, solver, Tc.geo->vertexDualAreas, 1)
                .nearestVertex();
        Tc.geo->unrequireVertexDualAreas();
    }

    VertexData<double> distancesToHorocycle =
        Tc.distanceOfHorocyclesTo(vInfinity);

    if (viz) {
        auto psMesh = polyscope::getSurfaceMesh("input_mesh");
        psMesh->addVertexCountQuantity(
            "vInfinity", std::vector<std::pair<size_t, int>>{
                             std::make_pair(vInfinity.getIndex(), 1)});
        psMesh->addVertexScalarQuantity("hyperbolic distance",
                                        distancesToHorocycle);
    }

    // ==========================================
    //     Compute uniformizing scale factors
    // ==========================================
    PetscWrapper optimizer(Tc, distancesToHorocycle, vInfinity.getIndex(),
                           std::numeric_limits<double>::infinity());

    // Set up petsc options
    std::vector<char*> argv;
    argv.emplace_back(
        const_cast<char*>("-tao_monitor")); // important for some reason
    argv.emplace_back(const_cast<char*>("-tao_type"));
    argv.emplace_back(const_cast<char*>("bnls"));
    // argv.emplace_back(const_cast<char*>("-citations"));
    // argv.emplace_back(const_cast<char*>("petsc_citations.txt"));

    std::clock_t start = std::clock();
    double duration;

    bool converged = optimizer.uniformize(argv.size(), argv.data());

    if (!converged) {
        std::cout << "Error: optimizer failed to converge" << vendl;
    }

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

    // ==========================================
    //              Map to sphere
    // ==========================================
    // Lay out in plane
    VertexData<Vector2> planeCoords = layOutTriangulation(Tc, vInfinity);

    // Project to sphere
    VertexData<Vector3> sphericalCoords(*Tc.mesh);
    auto stereographicProject = [](Vector2 v) -> Vector3 {
        double denom = 1 + v.x * v.x + v.y * v.y;
        return Vector3{2 * v.x / denom, 2 * v.y / denom,
                       (-1 + v.x * v.x + v.y * v.y) / denom};
    };

    for (Vertex v : Tc.mesh->vertices()) {
        sphericalCoords[v] = stereographicProject(planeCoords[v]);
    }
    sphericalCoords[vInfinity] = Vector3{0, 0, 1};

    // Apply Mobius centering
    Tb.geo->requireVertexDualAreas();
    VertexData<double> dualAreas =
        Tb.geo->vertexDualAreas.reinterpretTo(*Tc.mesh);
    sphericalCoords = mobiusCenter(*Tc.mesh, sphericalCoords, dualAreas);
    Tb.geo->unrequireVertexDualAreas();

    // ==========================================
    //             Common Refinement
    // ==========================================

    // Compute log scale factors for map from Tb -> sphere

    // sphere-inscribed edge lengths
    EdgeData<double> dstLen(*Tc.mesh);
    for (Edge e : Tc.mesh->edges()) {
        dstLen[e] = (sphericalCoords[e.firstVertex()] -
                     sphericalCoords[e.secondVertex()])
                        .norm();
    }

    VertexData<double> empiricalScaleFactors(*Tc.mesh, 0);
    const EdgeData<double>& srcLen = Tc.originalEdgeLengths;
    for (Vertex v : Tc.mesh->vertices()) {
        Halfedge he = v.halfedge();

        double srcL1   = srcLen[he.edge()];
        double srcLopp = srcLen[he.next().edge()];
        double srcL2   = srcLen[he.next().next().edge()];
        double dstL1   = dstLen[he.edge()];
        double dstLopp = dstLen[he.next().edge()];
        double dstL2   = dstLen[he.next().next().edge()];

        double expU = (dstL1 / srcL1) * (dstL2 / srcL2) / (dstLopp / srcLopp);

        empiricalScaleFactors[v] = log(expU);
    }
    Tc.setScaleFactors(empiricalScaleFactors);

    SphericalUniformizationResult result;
    std::tie(result.mesh, result.param, result.parentMap,
             result.interpolationMatrix) =
        computeCommonRefinementAndMatrix(Ta, Tb, Tc, vertexPositions,
                                         sphericalCoords, frontFaces);

    if (viz) {
        auto psSphere = polyscope::registerSurfaceMesh(
            "sphereMesh", sphericalCoords, Tc.mesh->getFaceVertexList());
        psSphere->setEnabled(false);

        auto psCommonRefinement = polyscope::registerSurfaceMesh(
            "CommonRefinement", result.mesh.vertexCoordinates,
            result.mesh.polygons);
        auto q = polyscope::addProjectiveSphericalParameterizationQuantity(
            *psCommonRefinement, "param", result.param);
        q->setEnabled(true);
    }

    return result;
}


static bool printedFNaN = false;
double f(PartiallyDecoratedTriangulation& tri, Face face) {
    double out = 0;
    tri.geo->requireCornerAngles();
    EdgeData<double>& lengths  = tri.geo->inputEdgeLengths;
    CornerData<double>& angles = tri.geo->cornerAngles;
    for (Halfedge he : face.adjacentHalfedges()) {
        double angle   = angles[he.next().next().corner()];
        double contrib = 0;
        contrib += angle * log(lengths[he.edge()]);
        contrib += lobachevsky(angle);

        if (std::isnan(contrib)) {
            if (!printedFNaN) {
                cerr << "f(Face) is NaN" << endl;
                cerr << "\t             angle: " << angle << endl;
                cerr << "\tlobachevsky(angle): " << lobachevsky(angle) << endl;
                cerr << "\t       log(length): " << log(lengths[he.edge()])
                     << endl;
                cerr << "\t edge lengths: ";
                for (Edge e : face.adjacentEdges()) {
                    cerr << lengths[e] << " ";
                }
                cerr << endl;
                printedFNaN = true;
            }
            return std::numeric_limits<double>::infinity();
            // my_assert(false, "objective error");
        }

        out += contrib;
    }
    if (out != out) out = std::numeric_limits<double>::infinity();
    return out;
}

double sphericalObjective(PartiallyDecoratedTriangulation& tri) {

    EdgeData<double>& lengths = tri.geo->inputEdgeLengths;

    // Vertex vInfinity = tri.mesh->vertex(655);
    // for (Vertex w : vInfinity.adjacentVertices()) {
    //     WATCH(w);
    //     std::cout << "\t" << w << ":\t" << 2 * PI * tri.logScaleFactors[w]
    //               << vendl;
    //     for (Edge e : w.adjacentEdges()) {
    //         if (tri.finite(e)) {
    //             std::cout << "\t" << src(e) << "-" << dst(e) << ":\t"
    //                       << 2 * PI * log(lengths[e]) << vendl;
    //         }
    //     }
    //     for (Face face : w.adjacentFaces()) {
    //         if (tri.finite(face)) {
    //             std::cout << "\t" << face.halfedge().vertex() << "-"
    //                       << face.halfedge().next().vertex() << "-"
    //                       << face.halfedge().next().next().vertex() << ":\t"
    //                       << 2 * f(tri, face) << vendl;
    //         }
    //     }
    // }

    double objective = 0;
    size_t iF        = 0;
    for (Face face : tri.mesh->faces()) {
        if (tri.finite(face)) {
            objective += 2 * f(tri, face);
            // if (iF++ < 5)
            //     std::cout << "f" << iF << ":\t" << 2 * f(tri, face) << vendl;
        }
    }

    size_t iE = 0;
    for (Edge e : tri.mesh->edges()) {
        if (tri.finite(e)) {
            objective -= PI * 2 * log(lengths[e]);
            // if (iE++ < 5)
            //     std::cout << "e" << iE << ":\t" << PI * 2 * log(lengths[e])
            //               << vendl;
        }
    }

    size_t iV = 0;
    for (Vertex v : tri.mesh->vertices()) {
        if (tri.finite(v)) {
            objective += 2 * PI * tri.logScaleFactors[v];
            // if (iV++ < 5)
            //     std::cout << "v" << iV << ":\t"
            //               << 2 * PI * tri.logScaleFactors[v] << vendl;
        }
    }

    // polyscope::getSurfaceMesh("input_mesh")
    //     ->addVertexScalarQuantity("u", tri.logScaleFactors);
    // polyscope::show();
    // throw_verbose_runtime_error("nyah");

    return objective;
}

VertexData<double>
sphericalGradientVertexData(PartiallyDecoratedTriangulation& tri) {
    VertexData<double> grad(*tri.mesh);
    VertexData<double> angleSums = tri.angleSums();
    for (Vertex v : tri.mesh->vertices()) {
        if (tri.finite(v)) {
            grad[v] = -angleSums[v] + PI * (tri.degF(v) - tri.degE(v) + 2);
        } else {
            grad[v] = 0;
        }
    }
    return grad;
}

Vector<double> sphericalGradient(PartiallyDecoratedTriangulation& tri) {
    return sphericalGradientVertexData(tri).toVector();
}

SparseMatrix<double> sphericalHessian(PartiallyDecoratedTriangulation& tri) {
    return tri.partialCotanLaplacian();
}

namespace ImplementationDetails {
VertexData<Vector2> layOutTriangulation(PartiallyDecoratedTriangulation& tri,
                                        Vertex vInfinity) {

    // Build Laplacian for mesh with vInfinity removed
    size_t nV               = tri.mesh->nVertices();
    VertexData<size_t> vIdx = tri.mesh->getVertexIndices();

    tri.geo->requireHalfedgeCotanWeights();
    std::vector<Eigen::Triplet<std::complex<double>>> TL;
    TL.emplace_back(vIdx[vInfinity], vIdx[vInfinity], 1);
    for (Halfedge he : tri.mesh->interiorHalfedges()) {
        if (!tri.finite(he.face())) continue;

        size_t iHead = vIdx[he.tipVertex()];
        size_t iTail = vIdx[he.tailVertex()];

        std::complex<double> weight = tri.geo->halfedgeCotanWeights[he];

        TL.emplace_back(iTail, iTail, weight);
        TL.emplace_back(iHead, iHead, weight);
        TL.emplace_back(iTail, iHead, -weight);
        TL.emplace_back(iHead, iTail, -weight);
    }
    for (size_t iC = 0; iC < nV; ++iC) TL.emplace_back(iC, iC, 1e-12);
    SparseMatrix<std::complex<double>> L(nV, nV);
    L.setFromTriplets(std::begin(TL), std::end(TL));

    // Build the mass matrix
    std::vector<Eigen::Triplet<std::complex<double>>> TM;
    tri.geo->requireFaceAreas();
    TM.emplace_back(vIdx[vInfinity], vIdx[vInfinity], 1);
    for (Face f : tri.mesh->faces()) {
        if (!tri.finite(f)) continue;
        std::complex<double> weight = tri.geo->faceAreas[f] / 3.;
        for (Vertex v : f.adjacentVertices()) {
            TM.emplace_back(vIdx[v], vIdx[v], weight / 3.);
        }
    }
    SparseMatrix<std::complex<double>> M(nV, nV);
    M.setFromTriplets(std::begin(TM), std::end(TM));

    // Build the area term
    std::complex<double> i(0, 1);
    std::vector<Eigen::Triplet<std::complex<double>>> TA;

    for (Halfedge heSpoke : vInfinity.outgoingHalfedges()) {
        Halfedge he = heSpoke.next().twin();
        size_t j    = vIdx[he.tipVertex()];
        size_t k    = vIdx[he.tailVertex()];

        TA.emplace_back(j, k, i * 0.25);
        TA.emplace_back(k, j, i * -0.25);
    }
    SparseMatrix<std::complex<double>> A(nV, nV);
    A.setFromTriplets(std::begin(TA), std::end(TA));

    SparseMatrix<std::complex<double>> EC = 0.5 * L - A;

    Vector<std::complex<double>> pos = ImplementationDetails::eig(EC, M, 1e-8);

    VertexData<Vector2> positions(*tri.mesh);
    for (Vertex v : tri.mesh->vertices()) {
        Vector2 cPos = Vector2::fromComplex(pos(vIdx[v]));
        positions[v] = cPos;
    }

    return positions;
}

Eigen::Vector3d toEigen(const Vector3& v) {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
}

Vector3 fromEigen(const Eigen::Vector3d& v) {
    return Vector3{v[0], v[1], v[2]};
}

// Returns vw^T
Eigen::Matrix3d outer(const Vector3& v, const Vector3& w) {
    return toEigen(v) * toEigen(w).transpose();
}

VertexData<Vector3> mobiusCenter(ManifoldSurfaceMesh& mesh,
                                 const VertexData<Vector3>& positions,
                                 VertexData<double> area, bool verbose) {

    double surfaceArea = area.toVector().sum();
    area /= surfaceArea;

    Vector3 centerOfMass;
    VertexData<Vector3> centeredPositions = positions;

    auto evalStep = [&](Vector3 c) {
        Vector3 newCenter{0, 0, 0};
        for (Vertex v : mesh.vertices()) {
            Vector3 p    = centeredPositions[v];
            Vector3 newP = (1 - norm2(c)) * (p + c) / norm2(p + c) + c;
            newCenter += area[v] * newP;
        }
        return norm(newCenter);
    };
    auto takeStep = [&](Vector3 c) {
        Vector3 newCenter{0, 0, 0};
        for (Vertex v : mesh.vertices()) {
            Vector3 p            = centeredPositions[v];
            centeredPositions[v] = (1 - norm2(c)) * (p + c) / norm2(p + c) + c;
            newCenter += area[v] * centeredPositions[v];
        }
        return newCenter;
    };

    centerOfMass    = takeStep(Vector3{0, 0, 0});
    double stepSize = 1;
    while (norm(centerOfMass) > 1e-3) {
        Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        for (Vertex v : mesh.vertices()) {
            Vector3 p = centeredPositions[v];
            J += 2 * area[v] * (I - outer(p, p));
        }

        Vector3 step = fromEigen(-J.inverse() * toEigen(centerOfMass));

        double oldObjective = norm(centerOfMass);


        // https://en.wikipedia.org/wiki/Golden-section_search
        double invphi  = 0.5 * (sqrt(5.) - 1); // 1/phi
        double invphi2 = 0.5 * (3 - sqrt(5.)); // 1/phi^2
        double a       = 0;
        double b       = 1;

        double h  = b - a;
        double c  = a + invphi2 * h;
        double d  = a + invphi * h;
        double fc = evalStep(c * step);
        double fd = evalStep(d * step);

        while (abs(c - d) > 1e-8) {
            if (fc < fd) {
                b  = d;
                d  = c;
                fd = fc;
                h  = invphi * h;
                c  = a + invphi2 * h;
                fc = evalStep(c * step);
            } else {
                a  = c;
                c  = d;
                fc = fd;
                h  = invphi * h;
                d  = a + invphi * h;
                fd = evalStep(d * step);
            }
        }
        double stepSize;
        if (fc < fd) {
            stepSize = (a + d) / 2;
        } else {
            stepSize = (c + b) / 2;
        }

        centerOfMass = takeStep(stepSize * step);

        if (verbose) {
            cout << "mu: " << centerOfMass << "\tnorm: " << norm(centerOfMass)
                 << "\tstepSize: " << stepSize << endl;
        }
    }
    return centeredPositions;
}

std::tuple<SimplePolygonMesh, std::vector<std::array<double, 4>>,
           VertexData<int>, SparseMatrix<double>>
computeCommonRefinementAndMatrix(Triangulation& Ta, Triangulation& Tb,
                                 Triangulation& Tc,
                                 const VertexData<Vector3>& initialPositions,
                                 const VertexData<Vector3>& sphericalPositions,
                                 const std::set<Face>& frontFaces) {

    using ImplementationDetails::computeCommonRefinement;
    using ImplementationDetails::Doublet;
    using ImplementationDetails::SparseVector;

    CornerData<DirectSum<Vector3, double>> homogSphericalPositions(*Tc.mesh);
    for (Corner c : Tc.mesh->corners()) {
        homogSphericalPositions[c] =
            DirectSum<Vector3, double>(sphericalPositions[c.vertex()], 1);
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
    std::vector<DirectSum<Vector3, double>> refinedHomogSphericalPositions;
    VertexData<int> refinedVertices;
    std::tie(commonRefinement, refinedInterpolation,
             refinedHomogSphericalPositions, refinedVertices) =
        computeCommonRefinement(Ta, Tb, Tc, trivialInterpolation,
                                homogSphericalPositions, isFrontFace);

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
            pos += d.val * initialPositions[Ta.mesh->vertex(d.col)];
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

    // Convert from CornerData to VertexData, and from DirectSum to std::array
    std::vector<std::array<double, 4>> niceRefinedSphereCoords(nV);
    size_t iC = 0;
    for (const std::vector<size_t>& face : commonRefinement.polygons) {
        for (size_t iV : face) {
            const DirectSum<Vector3, double>& q =
                refinedHomogSphericalPositions[iC];
            niceRefinedSphereCoords[iV] = {-q.x.x, q.x.y, q.x.z, q.y};
            iC++;
        }
    }

    return std::tie(commonRefinement, niceRefinedSphereCoords, refinedVertices,
                    interpolationMatrix);
}

} // namespace ImplementationDetails

} // namespace CEPS
