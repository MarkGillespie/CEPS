#pragma once
#include "Triangulation.h"
#include "Utils.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {

//== Optimizer for performing Euclidean uniformization/cone flattening
class Optimizer {
  public:
    Optimizer(Triangulation& tri, const VertexData<double>& targetAngleDefects);

    VertexData<double> targetAngleDefects;

    SparseMatrix<double> hessian();
    Vector<double> gradient(std::vector<Vertex> fixed);
    double objective();

    bool takeStep(Vector<double> stepDir, Vector<double> grad);

    bool uniformize(std::vector<Vertex> fixed, double tol = 1e-5);

    bool verbose      = false;
    bool extraVerbose = false;

  protected:
    Triangulation& tri;
    double f(Face f);
    VertexData<size_t> vIdx;
};
} // namespace CEPS
