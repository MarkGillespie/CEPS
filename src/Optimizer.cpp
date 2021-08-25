#include "Optimizer.h"

namespace CEPS {
Optimizer::Optimizer(Triangulation& tri_,
                     const VertexData<double>& targetAngleDefects_)
    : tri(tri_) {
    targetAngleDefects = targetAngleDefects_.reinterpretTo(*tri.mesh);
    vIdx               = tri.mesh->getVertexIndices();
}

SparseMatrix<double> Optimizer::hessian() {
    tri.geo->requireCotanLaplacian();
    SparseMatrix<double> H = tri.geo->cotanLaplacian;
    tri.geo->unrequireCotanLaplacian();
    return H;
}

Vector<double> Optimizer::gradient(std::vector<Vertex> fixed) {
    tri.geo->requireVertexGaussianCurvatures();
    const VertexData<double>& curvature = tri.geo->vertexGaussianCurvatures;
    Vector<double> grad(tri.mesh->nVertices());
    for (Vertex v : tri.mesh->vertices()) {
        grad[vIdx[v]] = curvature[v] - targetAngleDefects[v];
    }
    for (Vertex v : fixed) {
        grad[vIdx[v]] = 0;
    }
    tri.geo->unrequireVertexGaussianCurvatures();
    return grad;
}


double Optimizer::objective() {
    double obj = 0;
    // Vertex term
    for (Vertex v : tri.mesh->vertices()) {
        obj += (2 * M_PI - targetAngleDefects[v]) * tri.logScaleFactors[v];
    }
    // Edge term
    const EdgeData<double>& length = tri.geo->inputEdgeLengths;
    for (Edge e : tri.mesh->edges()) {
        obj -= 2 * M_PI * log(length[e]);
    }
    // Face term
    for (Face face : tri.mesh->faces()) {
        obj += 2 * f(face);
    }
    return obj;
}

bool Optimizer::takeStep(Vector<double> stepDirVec, Vector<double> gradVec) {
    VertexData<double> stepDir(*tri.mesh, stepDirVec);
    VertexData<double> oldU = tri.logScaleFactors;
    double oldF             = objective();
    double targetDecrease   = stepDirVec.dot(gradVec);

    if (verbose) {
        std::cout << "\t" << targetDecrease << "\t"
                  << gradVec.lpNorm<Eigen::Infinity>()
                  << "\t l2: " << gradVec.lpNorm<2>();
    }

    double stepInfinityNorm =
        fmax(stepDirVec.maxCoeff(), -stepDirVec.minCoeff());

    double stepSize = fmin(8, 30. / stepInfinityNorm);
    double newF     = oldF;

    size_t iter = 0;
    while (newF > oldF + 1e-2 * stepSize * targetDecrease && iter < 50) {
        stepSize /= 2;
        tri.setScaleFactors(oldU + stepSize * stepDir);
        tri.flipToDelaunay(GeometryType::HYPERBOLIC);
        newF = objective();
        iter += 1;
    }

    if (verbose) {
        std::cout << "\t" << stepSize
                  << "\tstepNorm (l_infty): " << (stepDirVec * stepSize).norm()
                  << std::endl;
    }

    return newF <= oldF + 1e-4 * stepSize * targetDecrease;
}

bool Optimizer::uniformize(std::vector<Vertex> fixed, double tol) {

    // Remove any duplicates
    std::sort(std::begin(fixed), std::end(fixed));
    fixed.erase(std::unique(std::begin(fixed), std::end(fixed)),
                std::end(fixed));

    tri.flipToDelaunay(GeometryType::HYPERBOLIC);
    Vector<double> grad = gradient(fixed);

    size_t flips = 0;

    // Selector for free vertices. Only used if some scale factors are
    // fixed
    Vector<bool> isFree;
    Vector<double> zeros;

    bool allFree = fixed.empty();

    Vector<double> ones;
    if (allFree) {
        ones = Vector<double>::Constant(tri.mesh->nVertices(),
                                        1. / tri.mesh->nVertices());
    } else {
        isFree = Vector<bool>::Ones(tri.mesh->nVertices(), true);
        zeros  = Vector<double>::Zero(fixed.size());
        for (Vertex v : fixed) isFree[vIdx[v]] = false;
    }


    if (verbose) {
        std::cout << "Initial objective: " << objective() << std::endl;
    }

    bool done         = grad.dot(grad) < tol;
    size_t totalSteps = 0;
    while (!done && totalSteps < 1000) {
        SparseMatrix<double> H = hessian();

        if (verbose) {
            std::cout << totalSteps << ": " << objective();
        }

        Vector<double> step;

        if (allFree) {
            step = solvePositiveDefinite(H, grad);

            // Project out constant part
            step -= step.sum() * ones;

            verbose_assert((H * step - grad).norm() < 1e-5, "solve failed");
            verbose_assert(abs(step.sum()) < 1e-5,
                           "step direction has constant part ?!");
            step = -step;
        } else {
            // If some scale factors are fixed, we need to project them
            // out of the solution
            BlockDecompositionResult<double> decomp =
                blockDecomposeSquare(H, isFree, false);
            Vector<double> rhsValsA, rhsValsB;
            decomposeVector(decomp, grad, rhsValsA, rhsValsB);
            Vector<double> freeResult =
                solvePositiveDefinite(decomp.AA, rhsValsA);
            step = reassembleVector(decomp, freeResult, zeros);
            step = -step;
        }

        done = (abs(step.dot(grad)) < tol);

        if (!takeStep(step, grad) && verbose) {
            std::cout << totalSteps << ". Line search failed\t" << objective()
                      << std::endl;
        }
        totalSteps++;

        // Update gradient
        grad = gradient(fixed);
    }

    if (verbose) {
        std::cout << "Finished in " << totalSteps << " steps." << std::endl;
        if (done) {
            std::cout << "Converged!" << std::endl;
        } else {
            std::cout << "Did not converge :'(" << std::endl;
        }
    }

    return done;
}

double Optimizer::f(Face face) {
    double out = 0;
    tri.geo->requireCornerAngles();
    EdgeData<double>& lengths  = tri.geo->inputEdgeLengths;
    CornerData<double>& angles = tri.geo->cornerAngles;
    for (Halfedge he : face.adjacentHalfedges()) {
        double angle   = angles[he.next().next().corner()];
        double contrib = 0;
        out += angle * log(lengths[he.edge()]);
        out += lobachevsky(angle);
    }
    if (out != out) out = std::numeric_limits<double>::infinity();
    tri.geo->unrequireCornerAngles();
    return out;
}
} // namespace CEPS
