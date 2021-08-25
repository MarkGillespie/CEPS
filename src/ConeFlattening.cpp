#include "ConeFlattening.h"

namespace CEPS {
ParameterizationResult parameterizeWithGreedyCones(ManifoldSurfaceMesh& mesh,
                                                   VertexPositionGeometry& geo,
                                                   bool viz,
                                                   bool checkInjectivity,
                                                   double uTol) {
    return parameterizeWithGivenCones(mesh, geo, placeCones(mesh, geo, uTol),
                                      viz, checkInjectivity);
}

ParameterizationResult parameterizeWithGivenCones(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
    const std::vector<Vertex>& coneVertices, bool viz, bool checkInjectivity) {
    std::vector<std::pair<Vertex, double>> prescribedScaleFactors;
    for (Vertex v : coneVertices) {
        prescribedScaleFactors.emplace_back(v, 0);
    }

    // Minimal area distortion boundary conditions.
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            prescribedScaleFactors.emplace_back(v, 0);
        }
    }

    return parameterize(mesh, geo, prescribedScaleFactors, {}, {}, viz,
                        checkInjectivity);
}

ParameterizationResult
parameterize(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
             std::vector<std::pair<Vertex, double>> prescribedScaleFactors,
             std::vector<std::pair<Vertex, double>> prescribedCurvatures,
             std::vector<size_t> imaginaryFaceIndices, bool viz,
             bool checkInjectivity) {
    if (mesh.nBoundaryLoops() == 0) { //== No boundary

        // == Identify front faces (all but the imaginary faces)
        std::set<size_t> frontFaceIndices;
        for (size_t iF = 0; iF < mesh.nFaces(); ++iF)
            frontFaceIndices.insert(iF);

        // Remove imaginary faces from front face list
        for (size_t iF : imaginaryFaceIndices) frontFaceIndices.erase(iF);

        return parameterizeHelper(mesh, geo, prescribedScaleFactors,
                                  prescribedCurvatures, frontFaceIndices, viz,
                                  checkInjectivity);

    } else { //== Mesh has boundary. We need to double the mesh

        //== Construct doubled mesh
        std::unique_ptr<ManifoldSurfaceMesh> doubledMesh;
        std::unique_ptr<VertexPositionGeometry> doubledGeometry;
        std::vector<Vertex> boundaryVs;
        VertexData<Vertex> twinMap;
        std::tie(doubledMesh, doubledGeometry, twinMap, boundaryVs) =
            doubleMesh(mesh, geo);

        //== Compute constraints on doubled mesh
        size_t nC = prescribedCurvatures.size();
        for (size_t iC = 0; iC < nC; ++iC) {
            std::pair<Vertex, double>& cone = prescribedCurvatures[iC];
            if (twinMap[cone.first] == cone.first) { // boundary vertices
                cone.second *= 2;
            } else { // interior vertices
                prescribedCurvatures.emplace_back(twinMap[cone.first],
                                                  cone.second);
            }
        }

        size_t nU = prescribedScaleFactors.size();
        for (size_t iC = 0; iC < nU; ++iC) {
            const std::pair<Vertex, double>& cone = prescribedScaleFactors[iC];
            if (twinMap[cone.first] == cone.first) { // boundary vertices
                                                     // (no action required)
            } else {                                 // interior vertices
                prescribedScaleFactors.emplace_back(twinMap[cone.first],
                                                    cone.second);
            }
        }

        // == Identify front faces
        std::set<size_t> frontFaceIndices;
        // By construction, the front faces are just the first half
        // of the list of faces

        for (size_t iF = 0; iF < mesh.nFaces(); ++iF)
            frontFaceIndices.insert(iF);

        // Remove imaginary faces from front face list
        for (size_t iF : imaginaryFaceIndices) frontFaceIndices.erase(iF);

        return parameterizeHelper(*doubledMesh, *doubledGeometry,
                                  prescribedScaleFactors, prescribedCurvatures,
                                  frontFaceIndices, viz, checkInjectivity);
    }
}

ParameterizationResult parameterizeHelper(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
    std::vector<std::pair<Vertex, double>> prescribedScaleFactors,
    std::vector<std::pair<Vertex, double>> prescribedCurvatures,
    std::set<size_t> frontFaceIndices, bool viz, bool checkInjectivity) {

    verbose_assert(mesh.nBoundaryLoops() == 0,
                   "parameterizeHelper only accepts closed meshes");

    // Check prescribed angle feasibility
    // These tests are not exhaustive - in actual meshes, angle defects cannot
    // be arbitrarily negative
    double prescribedCurvatureSum = 0;
    for (const auto cone : prescribedCurvatures)
        prescribedCurvatureSum += cone.second;
    double trueCurvatureSum = 2 * M_PI * mesh.eulerCharacteristic();

    double missingAngleSum = trueCurvatureSum - prescribedCurvatureSum;
    verbose_assert(
        abs(missingAngleSum) < 1e-6 ||
            (missingAngleSum < 1e-6 && !prescribedScaleFactors.empty()) ||
            (missingAngleSum <=
             2 * M_PI * prescribedScaleFactors.size() + 1e-6),
        "invalid boundary conditions: Gauss-Bonnet cannot be satisfied");

    Triangulation Ta(mesh, geo);
    VertexData<Vector3> vertexPositions =
        geo.inputVertexPositions.reinterpretTo(*Ta.mesh);
    std::set<Face> frontFaces;
    for (size_t iF : frontFaceIndices) frontFaces.insert(Ta.mesh->face(iF));


    // Copy Ta into Tb and then flip to intrinsic Delaunay
    Triangulation Tb(Ta, false);
    Tb.flipToDelaunay(GeometryType::EUCLIDEAN);

    // Copy Tb into Tc
    Triangulation Tc(Tb, true);

    // Initialize target angle defect to 0 at all vertices
    VertexData<double> targetAngleDefects(*Tc.mesh, 0);

    // Set target angle defect at prescribed vertices
    for (const std::pair<Vertex, double>& c : prescribedCurvatures) {
        targetAngleDefects[Tc.mesh->vertex(c.first.getIndex())] += c.second;
    }

    // Identify vertices with prescribed scale factors on triangulation Tc and
    // apply the desired scale factors
    VertexData<double> u(*Tc.mesh, 0);
    std::vector<Vertex> fixedVertices;
    for (const std::pair<Vertex, double>& c : prescribedScaleFactors) {
        Vertex v = Tc.mesh->vertex(c.first.getIndex());
        u[v]     = c.second;
        fixedVertices.push_back(v);
    }
    Tc.setScaleFactors(u);

    Optimizer newton(Tc, targetAngleDefects);

    // Compute uniformizing scale factors, fixing u at the fixed vertices
    bool converged = newton.uniformize(fixedVertices);

    if (!converged) {
        std::cout << "Error: optimizer failed to converge" << vendl;
    }

    // Identify set of all cones on Tc
    std::vector<Vertex> conesTc = fixedVertices;
    for (const std::pair<Vertex, double>& cone : prescribedCurvatures) {
        Vertex v = Tc.mesh->vertex(cone.first.getIndex());
        conesTc.push_back(v);
    }

    // Lay out triangles in plane to compute texture coordinates.
    // Normalize layout to have surface area 1
    CornerData<Vector2> textureCoords = layOutTriangulation(Tc, conesTc, 1);


    // Assemble the output data
    ParameterizationResult result;

    // Check for flipped triangles
    if (checkInjectivity) {
        std::tie(result.nFlippedTriangles, result.nZeroAreaTriangles) =
            checkTriangleOrientations(*Tc.mesh, textureCoords);
    }

    // Compute solution on common refinement, storing data in result
    std::tie(result.mesh, result.param, result.parentMap,
             result.interpolationMatrix) =
        computeCommonRefinementAndMatrix(Ta, Tb, Tc, vertexPositions,
                                         textureCoords, frontFaces);

    // TODO: convert parentmap to a map on the input mesh
    // TODO: fix interpolation matrix to work on doubled/filled meshes

    // Visualize solution
    if (viz) {
        auto psCommonRefinement = polyscope::registerSurfaceMesh(
            "Common Refinement", result.mesh.vertexCoordinates,
            result.mesh.polygons);

        auto pTex = polyscope::addProjectiveParameterizationQuantity(
            *psCommonRefinement, "projective texture", result.param);
        pTex->setEnabled(true);

        Vector<double> interpolatedScaleFactors =
            result.interpolationMatrix * Tc.logScaleFactors.toVector();
        psCommonRefinement->addVertexScalarQuantity("log scale factors",
                                                    interpolatedScaleFactors);

        // Measure cone angles and re-index using common refinement vertex
        // indices
        Tc.geo->requireVertexGaussianCurvatures();
        std::vector<std::pair<size_t, double>> refinedCones;
        VertexData<int> intrinsicParentMap =
            result.parentMap.reinterpretTo(*Tc.mesh);
        for (Vertex v : conesTc) {
            if (intrinsicParentMap[v] >= 0) {
                refinedCones.emplace_back(intrinsicParentMap[v],
                                          Tc.geo->vertexGaussianCurvatures[v]);
            }
        }
        Vector<double> interpolatedK =
            result.interpolationMatrix *
            Tc.geo->vertexGaussianCurvatures.toVector();
        psCommonRefinement->addVertexScalarQuantity("K", interpolatedK);
        Tc.geo->unrequireVertexGaussianCurvatures();

        auto coneQ = psCommonRefinement->addVertexIsolatedScalarQuantity(
            "cones", refinedCones);
        coneQ->setEnabled(true);

        // Ensure that the visualization range is symmetric, so that positive
        // cones are red and negative cones are blue when visualized with the
        // coolwarm colormap
        symmetrizeVizRange(*coneQ);
    }

    return result;
}

} // namespace CEPS
