#include "GeometryUtils.h"

namespace CEPS {

// stolen from geometrycentral/surface/direction_fields.h
VertexData<int> computeVertexIndex(ManifoldSurfaceMesh& mesh,
                                   IntrinsicGeometryInterface& geo,
                                   const FaceData<Vector2>& directionField,
                                   int nSym) {
    geo.requireTransportVectorsAcrossHalfedge();
    geo.requireVertexAngleSums();
    // Store the result here
    VertexData<int> indices(mesh);
    for (Vertex v : mesh.vertices()) {
        // Trace the direction field around the face and see how many times it
        // spins!
        double totalRot = 0;
        for (Halfedge he : v.incomingHalfedges()) {
            if (he.edge().isBoundary()) continue;

            // Compute the rotation along the halfedge implied by the field
            Vector2 x0 = directionField[he.face()].pow(nSym);
            Vector2 x1 = directionField[he.twin().face()].pow(nSym);
            Vector2 transport =
                geo.transportVectorsAcrossHalfedge[he].pow(nSym);
            // Find the difference in angle
            double theta0     = arg(transport * x0);
            double theta1     = arg(x1);
            double deltaTheta = arg(x1 / (transport * x0));
            totalRot += deltaTheta;
        }
        double angleDefect = (v.isBoundary())
                                 ? 1. * M_PI - geo.vertexAngleSums[v]
                                 : 2. * M_PI - geo.vertexAngleSums[v];
        totalRot += angleDefect * nSym;
        // Compute the net rotation and corresponding index
        // should be very close to a multiple of 2PI
        int index  = static_cast<int>(std::round(totalRot / (2 * PI)));
        indices[v] = index;
    }
    return indices;
}

std::vector<std::pair<Vertex, double>>
lumpCones(ManifoldSurfaceMesh& mesh,
          const std::vector<std::pair<Vertex, double>>& cones) {

    std::vector<Vertex> coneVs;
    VertexData<double> coneAngle(mesh, 0);
    for (const std::pair<Vertex, double>& cone : cones) {
        coneVs.push_back(cone.first);
        coneAngle[cone.first] = cone.second;
    }

    // Distribute angles from any cones with angle >= 2 pi
    // Note that this is not guaranteed to work with really bad cone
    // configurations. But it generally seems to work well enough
    for (Vertex v : coneVs) {
        while (coneAngle[v] >= 2 * M_PI) {
            Vertex smallestNeighbor;
            for (Vertex w : v.adjacentVertices()) {
                if (smallestNeighbor == Vertex() ||
                    coneAngle[w] < coneAngle[smallestNeighbor]) {
                    smallestNeighbor = w;
                }
            }
            if (coneAngle[smallestNeighbor] < -1e-8) {
                coneAngle[v] += coneAngle[smallestNeighbor];
                coneAngle[smallestNeighbor] = 0;
            } else if (abs(coneAngle[smallestNeighbor]) < 1e-8) {
                // If neighbor is not a cone yet, add it to the cone list
                coneAngle[v] -= M_PI / 4;
                coneAngle[smallestNeighbor] += M_PI / 4;
                coneVs.push_back(smallestNeighbor);
            } else {
                coneAngle[v] -= M_PI / 4;
                coneAngle[smallestNeighbor] += M_PI / 4;
            }
        }
    }

    // Merge neighboring cones
    bool done = false;
    while (!done) {
        done = true;
        for (Vertex cone : coneVs) {
            if (abs(coneAngle[cone]) < 1e-12) continue;

            for (Vertex w : cone.adjacentVertices()) {
                if (abs(coneAngle[w]) > 1e-12 &&
                    coneAngle[cone] + coneAngle[w] < 2 * M_PI * 0.8) {
                    coneAngle[cone] += coneAngle[w];
                    coneAngle[w] = 0;
                    done         = false;
                }
            }
        }
    }

    std::vector<std::pair<Vertex, double>> lumpedCones;
    for (Vertex v : coneVs) {
        if (abs(coneAngle[v]) > 1e-12) {
            lumpedCones.push_back(std::make_pair(v, coneAngle[v]));
        }
    }
    return lumpedCones;
}

// Doubles a mesh and geometry, gluing two copies of the mesh along their
// boundary to produce a single mesh without boundary.
// Returns the new mesh, the new geometry, the doubling map (taking each vertex
// to its twin on the other copy of the original mesh), and a list of boundary
// vertices. The doubling map sends boundary vertices to themselves
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
           std::unique_ptr<VertexPositionGeometry>, VertexData<Vertex>,
           std::vector<Vertex>>
doubleMesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {
    std::unique_ptr<ManifoldSurfaceMesh> doubledMesh;
    VertexData<Vertex> parentVtx;
    std::vector<Vertex> boundaryVertices;

    size_t nV = mesh.nVertices();
    size_t nB = nV - mesh.nInteriorVertices();

    // Double the mesh
    std::tie(doubledMesh, parentVtx, boundaryVertices) =
        ImplementationDetails::doubleMesh(mesh);

    // Copy over vertex positions to the doubled vertices
    VertexData<Vector3> doubledVertexPositions(*doubledMesh);
    for (Vertex v : doubledMesh->vertices()) {
        doubledVertexPositions[v] = geo.inputVertexPositions[parentVtx[v]];
    }

    // Create the doubled geometry
    VertexPositionGeometry* doubledGeo =
        new VertexPositionGeometry(*doubledMesh, doubledVertexPositions);

    // For each vertex in the original mesh, record its children in the doubled
    // mesh
    VertexData<std::vector<Vertex>> children(mesh);
    for (Vertex v : doubledMesh->vertices()) {
        children[parentVtx[v]].push_back(v);
    }

    // The twin map for doubled vertices just maps any original parent vertex's
    // children to each other
    VertexData<Vertex> twin(*doubledMesh);
    for (Vertex v : mesh.vertices()) {
        if (children[v].size() == 1) {
            twin[children[v][0]] = children[v][0];
        } else if (children[v].size() == 2) {
            twin[children[v][0]] = children[v][1];
            twin[children[v][1]] = children[v][0];
        } else if (children[v].size() == 0) {
            throw_verbose_runtime_error("doubling mesh lost a vertex");
        } else {
            throw_verbose_runtime_error(
                "doubling mesh produced more than two copies of a "
                "vertex.");
        }
    }

    return std::make_tuple(std::move(doubledMesh),
                           std::unique_ptr<VertexPositionGeometry>(doubledGeo),
                           twin, boundaryVertices);
}

std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {
    std::unique_ptr<ManifoldSurfaceMesh> meshCopy = mesh.copy();
    geo.requireEdgeLengths();
    EdgeData<double> lengthsCopy = geo.edgeLengths.reinterpretTo(*meshCopy);
    std::unique_ptr<EdgeLengthGeometry> geoCopy =
        std::make_unique<EdgeLengthGeometry>(*meshCopy, lengthsCopy);
    return std::make_pair(std::move(meshCopy), std::move(geoCopy));
}


std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, EdgeLengthGeometry& geo) {
    std::unique_ptr<ManifoldSurfaceMesh> meshCopy = mesh.copy();
    EdgeData<double> lengthsCopy =
        geo.inputEdgeLengths.reinterpretTo(*meshCopy);
    std::unique_ptr<EdgeLengthGeometry> geoCopy =
        std::make_unique<EdgeLengthGeometry>(*meshCopy, lengthsCopy);
    return std::make_pair(std::move(meshCopy), std::move(geoCopy));
}

std::pair<size_t, size_t> checkTriangleOrientations(ManifoldSurfaceMesh& mesh,
                                                    CornerData<Vector2>& uv) {
    size_t nFlipped  = 0;
    size_t nZeroArea = 0;

    ImplementationDetails::exactinit(); // initialize predicates.c
    for (Face f : mesh.faces()) {
        Vector2 p = uv[f.halfedge().corner()];
        Vector2 q = uv[f.halfedge().next().corner()];
        Vector2 r = uv[f.halfedge().next().next().corner()];

        double triOrientation = ImplementationDetails::orientation(p, q, r);
        if (triOrientation < 0) {
            nFlipped++;
        } else if (triOrientation == 0) {
            nZeroArea++;
        }
    }

    return std::make_pair(nFlipped, nZeroArea);
}

namespace ImplementationDetails {

// Doubles a mesh, gluing two copies of the mesh along their
// boundary to produce a single mesh without boundary.
// Returns the new mesh, the parent vertex in the original mesh for each doubled
// vertex, and a list of boundary vertices
// Beware - the gluing involves tricky halfedge manipulations
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, VertexData<Vertex>,
           std::vector<Vertex>>
doubleMesh(ManifoldSurfaceMesh& mesh) {
    // twin generation code stolen from surface_mesh.cpp in geometry-central


    size_t nV               = mesh.nVertices();
    size_t nF               = mesh.nFaces();
    VertexData<size_t> vIdx = mesh.getVertexIndices();

    // Mark boundary vertices
    std::vector<char> isBdy(nV, false);
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isBdy[vIdx[v]] = true;
        }
    }

    // Extract face list
    std::vector<std::vector<size_t>> frontFaces = mesh.getFaceVertexList();

    // Double faces (but leave boundary vertices glued)
    std::vector<std::vector<size_t>> backFaces;
    backFaces.reserve(nF);
    for (std::vector<size_t> face : frontFaces) {
        std::reverse(std::begin(face), std::end(face));
        for (size_t& iV : face) {
            if (!isBdy[iV]) {
                iV += nV;
            }
        }
        backFaces.push_back(face);
    }

    // =====================================================================
    //                 Compute the Doubled Mesh's Twin Maps
    // =====================================================================
    // Note that the twin maps on the back faces are different since the back
    // faces have their orientation flipped
    FaceData<size_t> fIdx = mesh.getFaceIndices();
    std::vector<std::vector<std::tuple<size_t, size_t>>> frontTwins(nF);
    std::vector<std::vector<std::tuple<size_t, size_t>>> backTwins(nF);

    HalfedgeData<size_t> iHeInFrontFace(mesh);
    for (Face f : mesh.faces()) {
        size_t iHe = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            iHeInFrontFace[he] = iHe;
            iHe++;
        }
    }

    for (Face f : mesh.faces()) {
        size_t iF = fIdx[f];
        size_t D  = f.degree();

        // Reverse orientation
        auto opp = [&](size_t i) { return (D + D - i - 2) % D; };

        std::vector<std::tuple<size_t, size_t>>& frontTwin = frontTwins[iF];
        std::vector<std::tuple<size_t, size_t>>& backTwin  = backTwins[iF];
        frontTwin.resize(D);
        backTwin.resize(D);

        size_t i = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            verbose_assert(i + 1 <= D, "face has too many halfedges");
            verbose_assert(i + 1 <= D && D < D + 1 + i, "I can't do algebra");
            if (he.edge().isBoundary()) {
                size_t heTIndFront = iHeInFrontFace[he];
                size_t heTIndBack  = opp(heTIndFront);
                frontTwin[i]       = std::make_tuple(iF + nF, heTIndBack);
                backTwin[opp(i)]   = std::make_tuple(iF, heTIndFront);
            } else {
                Halfedge heT       = he.twin();
                size_t fT          = fIdx[heT.face()];
                size_t heTIndFront = iHeInFrontFace[heT];
                size_t heTIndBack  = opp(heTIndFront);
                frontTwin[i]       = std::make_tuple(fT, heTIndFront);
                backTwin[opp(i)]   = std::make_tuple(fT + nF, heTIndBack);
            }
            i++;
        }
    }


    std::vector<std::vector<size_t>> faceList = frontFaces;
    faceList.insert(std::end(faceList), std::begin(backFaces),
                    std::end(backFaces));
    std::vector<size_t> newVInd = stripUnusedVertices(faceList, 2 * nV);

    std::vector<std::vector<std::tuple<size_t, size_t>>> twins = frontTwins;
    twins.insert(std::end(twins), std::begin(backTwins), std::end(backTwins));

    // Construct a mesh from these faces and twin maps
    ManifoldSurfaceMesh* doubledMesh = new ManifoldSurfaceMesh(faceList, twins);

    // Record the doubled mesh's boundary
    std::vector<Vertex> boundaryVertices;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            boundaryVertices.push_back(doubledMesh->vertex(vIdx[v]));
        }
    }

    // Compute a parent in the original mesh for each vertex of the doubled mesh
    VertexData<Vertex> parentVtx(*doubledMesh);
    for (Vertex v : mesh.vertices()) {
        // Original vertices are their own parent
        parentVtx[doubledMesh->vertex(vIdx[v])] = v;
        if (!v.isBoundary())
            parentVtx[doubledMesh->vertex(newVInd[vIdx[v] + nV])] = v;
    }

    return std::make_tuple(std::unique_ptr<ManifoldSurfaceMesh>(doubledMesh),
                           parentVtx, boundaryVertices);
}

// Take in a list of faces, and an upper bound on the number vertices. Reindex
// the faces to remove any unused vertices and return the map sending an old
// index to its new compressed value
// Stolen from geometrycentral/surface/simple_polygon_mesh.cpp
std::vector<size_t>
stripUnusedVertices(std::vector<std::vector<size_t>>& faceList, size_t nV) {
    // Check which indices are used
    std::vector<char> vertexUsed(nV, false);
    for (const std::vector<size_t>& face : faceList) {
        for (size_t iV : face) {
            verbose_assert(iV < nV, "invalid vertex index " +
                                        std::to_string(iV) +
                                        " (there should only be " +
                                        std::to_string(nV) + " vertices)");
            vertexUsed[iV] = true;
        }
    }


    // Re-index
    std::vector<size_t> newInd(nV, INVALID_IND);
    size_t nNewV = 0;
    for (size_t iOldV = 0; iOldV < nV; iOldV++) {
        if (!vertexUsed[iOldV]) continue;
        size_t iNewV  = nNewV++;
        newInd[iOldV] = iNewV;
    }

    // Translate the polygon listing
    for (std::vector<size_t>& face : faceList) {
        for (auto& iV : face) {
            iV = newInd[iV];
        }
    }

    return newInd;
}

double orientation(const Vector2& p, const Vector2& q, const Vector2& r) {
    // if (!predicatesInitialized) {
    //     exactinit();
    //     predicatesInitialized = true;
    // }
    std::array<double, 6> values{p.x, p.y, q.x, q.y, r.x, r.y};
    return orient2d(&values[0], &values[2], &values[4]);
}


// Computes the smallest eigenvector of M^-1*A orthogonal to 1
Vector<std::complex<double>> eig(SparseMatrix<std::complex<double>>& A,
                                 const SparseMatrix<std::complex<double>>& M,
                                 double tol) {

    Vector<std::complex<double>> ones =
        Vector<std::complex<double>>::Ones(A.rows());
    auto norm = [&](const Vector<std::complex<double>>& v) {
        return std::sqrt(std::abs(v.dot(M * v)));
    };
    ones /= norm(ones);

    size_t N = A.rows();
    PositiveDefiniteSolver<std::complex<double>> solver(A);

    auto projectOutOnes = [&](Vector<std::complex<double>>& x) {
        std::complex<double> proj = ((ones.dot(M * x)));
        x -= proj * ones;
    };
    auto residual = [&](const Vector<std::complex<double>>& v) {
        std::complex<double> candidateEigenvalue = v.dot(A * v);
        Vector<std::complex<double>> err = A * v - candidateEigenvalue * M * v;
        return norm(err) / norm(v);
    };

    Vector<std::complex<double>> u = Vector<std::complex<double>>::Random(N);
    projectOutOnes(u);
    Vector<std::complex<double>> x = u;
    size_t iter                    = 0;
    while (residual(x) > tol && iter++ < 1000) {
        // Solve
        solver.solve(x, M * u);

        projectOutOnes(x);
        x /= norm(x);

        // Update
        u = x;
    }

    return x;
}

} // namespace ImplementationDetails
} // namespace CEPS
