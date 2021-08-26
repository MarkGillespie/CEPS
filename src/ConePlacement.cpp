#include "ConePlacement.h"

namespace CEPS {
std::vector<Vertex> placeCones(ManifoldSurfaceMesh& mesh,
                               IntrinsicGeometryInterface& geo, double uTol,
                               size_t minNumCones, size_t gaussSeidelIterations,
                               bool verbose) {

    // We only place cones on interior vertices. If there aren't any, don't
    // place any cones
    if (mesh.nInteriorVertices() == 0) return {};

    // Using the ordinary cotan Laplacian leads to numerical issues on
    // near-degenerate models (e.g. thingi10k 424212_7.ply)
    SparseMatrix<double> L =
        ImplementationDetails::intrinsicDelaunayLaplacian(mesh, geo);


    // Also regularize slightly
    Eigen::SparseMatrix<double> eye(mesh.nVertices(), mesh.nVertices());
    eye.setIdentity();

    L = L + 1e-6 * eye;

    auto getWorstVertex = [&](const VertexData<double>& vec) -> Vertex {
        Vertex worstSoFar = mesh.vertex(0);
        for (Vertex v : mesh.vertices()) {
            if (abs(vec[v]) > abs(vec[worstSoFar])) worstSoFar = v;
        }
        return worstSoFar;
    };

    if (verbose) std::cout << "\t computing initial greedy cones" << std::endl;

    // Place first cone by approximately uniformizing to a constant curvature
    // metric and picking the vertex with worst scale factor (We would flatten,
    // except you can't necessarily flatten with one cone)
    VertexData<double> scaleFactors =
        ImplementationDetails::computeConstantCurvatureScaleFactors(mesh, geo,
                                                                    L);
    Vertex worstVertex = getWorstVertex(scaleFactors);

    std::vector<Vertex> cones;
    cones.push_back(worstVertex);

    // Place cones by approximately flattening and picking the vertex
    // with worst scale factor. Keep doing this until you have at least
    // minNumCones cones, and the worst scale factor is at most uTol
    while (cones.size() < minNumCones ||
           abs(scaleFactors[worstVertex]) > uTol) {
        scaleFactors =
            ImplementationDetails::computeScaleFactors(mesh, geo, L, cones);
        worstVertex = getWorstVertex(scaleFactors);
        cones.push_back(worstVertex);
    }

    if (verbose) std::cout << "\t improving cones" << std::endl;

    // Improve cones by removing each cone and replacing it by the vertex with
    // the biggest scale factor
    for (size_t iGS = 0; iGS < gaussSeidelIterations; ++iGS) {
        for (size_t iCone = 0; iCone < cones.size(); ++iCone) {
            VertexData<double> scaleFactors =
                ImplementationDetails::computeScaleFactors(mesh, geo, L, cones,
                                                           iCone);
            cones[iCone] = getWorstVertex(scaleFactors);
        }
    }

    return cones;
}

namespace ImplementationDetails {

// Compute scale factors which give the mesh constant curvature
VertexData<double>
computeConstantCurvatureScaleFactors(ManifoldSurfaceMesh& mesh,
                                     IntrinsicGeometryInterface& geo,
                                     const SparseMatrix<double>& L) {
    geo.requireCotanLaplacian();
    geo.requireVertexGaussianCurvatures();

    Vector<double> weightedCurvature = -geo.vertexGaussianCurvatures.toVector();

    double totalCurvature = weightedCurvature.sum();
    double avgCurvature   = totalCurvature / (double)mesh.nVertices();

    Vector<double> constantCurvature =
        Vector<double>::Constant(mesh.nVertices(), avgCurvature);

    Vector<double> curvatureDifference = weightedCurvature - constantCurvature;


    // Identify interior vertices
    VertexData<size_t> vIdx = mesh.getVertexIndices();
    Vector<bool> isInterior = Vector<bool>::Constant(mesh.nVertices(), true);

    size_t nB = 0;
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior[vIdx[v]] = false;
            nB++;
        }
    }


    // Solve Lu = curvatureDifference w/ zero Dirichlet boundary conditions
    BlockDecompositionResult<double> laplacianBlocks =
        blockDecomposeSquare(L, isInterior, true);

    Vector<double> boundaryCurvature, interiorCurvature;
    decomposeVector(laplacianBlocks, weightedCurvature, interiorCurvature,
                    boundaryCurvature);

    Vector<double> interiorResult =
        solvePositiveDefinite(laplacianBlocks.AA, interiorCurvature);
    Vector<double> zero = Eigen::VectorXd::Zero(nB);
    Vector<double> totalResult =
        reassembleVector(laplacianBlocks, interiorResult, zero);

    return VertexData<double>(mesh, totalResult);
}

// Compute scale factors which flatten mesh using the given cones.
// If vSkip1 or vSkip2 are nonzero, ignore the corresponding cones (i.e. flatten
// those vertices too)
VertexData<double> computeScaleFactors(ManifoldSurfaceMesh& mesh,
                                       IntrinsicGeometryInterface& geo,
                                       const SparseMatrix<double>& L,
                                       const std::vector<Vertex> cones,
                                       int vSkip1, int vSkip2) {
    geo.requireCotanLaplacian();
    geo.requireVertexGaussianCurvatures();

    // interior vertices are the ones that aren't cones
    VertexData<size_t> vIdx = mesh.getVertexIndices();
    Vector<bool> isInterior = Vector<bool>::Constant(mesh.nVertices(), true);

    size_t nCones = 0;
    for (size_t iV = 0; iV < cones.size(); ++iV) {
        if ((int)iV != vSkip1 && (int)iV != vSkip2) {
            if (isInterior[vIdx[cones[iV]]]) nCones++;
            isInterior[vIdx[cones[iV]]] = false;
        }
    }
    // If the mesh has boundary, make every boundary vertex a cone
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex v : b.adjacentVertices()) {
            isInterior[vIdx[v]] = false;
            nCones++;
        }
    }

    BlockDecompositionResult<double> laplacianBlocks =
        blockDecomposeSquare(L, isInterior, true);

    Vector<double> weightedCurvature = -geo.vertexGaussianCurvatures.toVector();

    Vector<double> coneCurvature, interiorCurvature;
    decomposeVector(laplacianBlocks, weightedCurvature, interiorCurvature,
                    coneCurvature);

    Vector<double> interiorResult =
        solvePositiveDefinite(laplacianBlocks.AA, interiorCurvature);
    Vector<double> zero = Eigen::VectorXd::Zero(nCones);
    Vector<double> totalResult =
        reassembleVector(laplacianBlocks, interiorResult, zero);

    return VertexData<double>(mesh, totalResult);
}

SparseMatrix<double> intrinsicDelaunayLaplacian(ManifoldSurfaceMesh& mesh,
                                                IntrinsicGeometryInterface& geo,
                                                double mollify) {
    std::unique_ptr<ManifoldSurfaceMesh> meshCpy = mesh.copy();

    geo.requireEdgeLengths();
    EdgeData<double> eLen = geo.edgeLengths.reinterpretTo(*meshCpy);
    geo.unrequireEdgeLengths();

    EdgeLengthGeometry eGeo(*meshCpy, eLen);

    mollifyIntrinsic(*meshCpy, eGeo.inputEdgeLengths, mollify);
    flipToDelaunay(*meshCpy, eGeo.inputEdgeLengths);

    eGeo.requireHalfedgeCotanWeights();
    for (Halfedge he : eGeo.mesh.interiorHalfedges()) {
        if (!std::isfinite(eGeo.halfedgeCotanWeights[he])) {
            std::cout << "Bad cotan weight: " << eGeo.halfedgeCotanWeights[he]
                      << vendl;
            std::array<double, 3> edgeLengths{
                eGeo.inputEdgeLengths[he.edge()],
                eGeo.inputEdgeLengths[he.next().edge()],
                eGeo.inputEdgeLengths[he.next().next().edge()]};
            std::cout << "\tedge lengths:";
            for (double l : edgeLengths) {
                std::cout << " " << l;
            }
            std::cout << vendl;
        }
    }


    eGeo.requireCotanLaplacian();

    return eGeo.cotanLaplacian;
}

} // namespace ImplementationDetails
} // namespace CEPS
