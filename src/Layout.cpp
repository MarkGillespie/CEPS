#include "Layout.h"

namespace CEPS {

CornerData<Vector2> layOutTriangulation(Triangulation& tri,
                                        const std::vector<Vertex>& cones,
                                        double targetSurfaceArea) {


    EdgeData<double> edgeCost(*tri.mesh, 1);
    for (Edge e : tri.mesh->edges()) verbose_assert(edgeCost[e] == 1, "??");
    return layOutTriangulation(tri, cones, edgeCost, targetSurfaceArea);
}

CornerData<Vector2> layOutTriangulation(Triangulation& tri,
                                        const std::vector<Vertex>& cones,
                                        const EdgeData<double>& edgeCutCosts,
                                        double targetSurfaceArea) {

    std::set<Edge> cut =
        ImplementationDetails::goodCuts(*tri.mesh, cones, edgeCutCosts);

    // Hack to visualize cuts for debugging
    // auto drawCut = [&](const std::set<Edge>& cut, std::string name) {
    //     polyscope::SurfaceMesh* psMesh =
    //         polyscope::getSurfaceMesh("input_mesh");

    //     std::vector<glm::vec3> cutVertexPositions;
    //     std::vector<std::array<size_t, 2>> cutEdges;
    //     for (Edge e : cut) {
    //         cutEdges.push_back(
    //             {cutVertexPositions.size(), cutVertexPositions.size() + 1});
    //         cutVertexPositions.push_back(
    //             psMesh->vertices[e.firstVertex().getIndex()]);
    //         cutVertexPositions.push_back(
    //             psMesh->vertices[e.secondVertex().getIndex()]);
    //     }
    //     polyscope::registerCurveNetwork(name, cutVertexPositions, cutEdges);
    // };
    // drawCut(cut, "cut");

    CornerData<size_t> cutCornerIndices;
    size_t nV;
    std::tie(cutCornerIndices, nV) =
        ImplementationDetails::indexCorners(*tri.mesh, cut);

    // Check that cut mesh is a topological disk
    int cutEulerCharacteristic = 0;
    cutEulerCharacteristic =
        nV + tri.mesh->nFaces() - tri.mesh->nEdges() - cut.size();
    WATCH(cutEulerCharacteristic);


    // Build Laplacian for the cut mesh
    // (We can't use geometry-central's Laplacian due to the cuts)
    tri.geo->requireHalfedgeCotanWeights();
    std::vector<Eigen::Triplet<std::complex<double>>> TL;
    for (Halfedge he : tri.mesh->interiorHalfedges()) {
        Corner cTail = he.corner();
        Corner cHead = he.next().corner();

        size_t iCHead = cutCornerIndices[cHead];
        size_t iCTail = cutCornerIndices[cTail];

        std::complex<double> weight = tri.geo->halfedgeCotanWeights[he];

        TL.emplace_back(iCTail, iCTail, weight);
        TL.emplace_back(iCHead, iCHead, weight);
        TL.emplace_back(iCTail, iCHead, -weight);
        TL.emplace_back(iCHead, iCTail, -weight);
    }
    for (size_t iC = 0; iC < nV; ++iC) TL.emplace_back(iC, iC, 1e-12);
    SparseMatrix<std::complex<double>> L(nV, nV);
    L.setFromTriplets(std::begin(TL), std::end(TL));

    // Build the mass matrix
    std::vector<Eigen::Triplet<std::complex<double>>> TM;
    tri.geo->requireFaceAreas();
    for (Face f : tri.mesh->faces()) {
        std::complex<double> weight = tri.geo->faceAreas[f] / 3.;
        for (Corner c : f.adjacentCorners()) {
            size_t iC = cutCornerIndices[c];
            TM.emplace_back(iC, iC, weight / 3.);
        }
    }
    SparseMatrix<std::complex<double>> M(nV, nV);
    M.setFromTriplets(std::begin(TM), std::end(TM));

    // Build the area term
    std::complex<double> i(0, 1);
    std::vector<Eigen::Triplet<std::complex<double>>> TA;

    for (Edge e : cut) {
        for (Halfedge he : e.adjacentInteriorHalfedges()) {
            size_t j = cutCornerIndices[he.next().corner()];
            size_t k = cutCornerIndices[he.corner()];

            TA.emplace_back(
                Eigen::Triplet<std::complex<double>>(j, k, i * 0.25));
            TA.emplace_back(
                Eigen::Triplet<std::complex<double>>(k, j, i * -0.25));
        }
    }
    SparseMatrix<std::complex<double>> A(nV, nV);
    A.setFromTriplets(std::begin(TA), std::end(TA));

    SparseMatrix<std::complex<double>> EC = 0.5 * L - A;

    Vector<std::complex<double>> pos = ImplementationDetails::eig(EC, M, 1e-8);

    CornerData<Vector2> positions(*tri.mesh);
    for (Vertex v : tri.mesh->vertices()) {
        for (Corner c : v.adjacentCorners()) {
            Vector2 cPos = Vector2::fromComplex(pos(cutCornerIndices[c]));
            positions[c] = cPos;
        }
    }

    // Normalize parameterization to have the desired surface area
    double layoutSurfaceArea = 0;
    for (Face f : tri.mesh->faces()) {
        Vector2 p = positions[f.halfedge().corner()];
        Vector2 q = positions[f.halfedge().next().corner()];
        Vector2 r = positions[f.halfedge().next().next().corner()];

        double fArea = abs(cross(q - p, r - p)) / 2.;
        layoutSurfaceArea += fArea;
    }

    positions *= sqrt(targetSurfaceArea / layoutSurfaceArea);

    return positions;
}
} // namespace CEPS
