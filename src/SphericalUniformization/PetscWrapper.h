#pragma once

#include <petsctao.h>

#include "PartiallyDecoratedTriangulation.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace CEPS;

using std::cerr;
using std::cout;
using std::endl;

namespace CEPS {

// Forward declarations of spherical energy, defined in SphericalUniformization
double sphericalObjective(PartiallyDecoratedTriangulation& tri);
Vector<double> sphericalGradient(PartiallyDecoratedTriangulation& tri);
SparseMatrix<double> sphericalHessian(PartiallyDecoratedTriangulation& tri);

class PetscWrapper {

  public:
    PetscWrapper(PartiallyDecoratedTriangulation& tri, size_t vInfinityIndex_);

    PetscWrapper(PartiallyDecoratedTriangulation& tri, size_t vInfinityIndex_,
                 double vInfinityScaling);

    PetscWrapper(PartiallyDecoratedTriangulation& tri,
                 VertexData<double>& distancesToHorocycle_,
                 size_t vInfinityIndex_, double vInfinityScaling);
    void setup(double vInfinityScaling);

    PartiallyDecoratedTriangulation& tri;
    Vertex vInfinity;
    size_t vInfinityIndex;
    VertexData<double> distancesToHorocycle;
    VertexData<size_t> vIdx;

    bool uniformize(int argc, char** argv);
    PetscErrorCode testOptimize(int argc, char** argv);

    void kktError(double& dualInfeasibility,
                  double& complementarySlacknessError, double& angleError);

    bool activeConstraint(Vertex v);
    void adjustTriangulation();
    void cleanUpTriangulation();

    bool useFlips   = true;
    double flipTime = 0;
};

/*
  FormFunctionGradient - Evaluates the function, f(X), and gradient, G(X).

  Input Parameters:
  .   tao  - the Tao context
  .   X    - input vector
  .   ptr  - optional user-defined context, as set by
  TaoSetFunctionGradient()

  Output Parameters:
  .   G - vector containing the newly evaluated gradient
  .   f - function value

  Note:
  Some optimization methods ask for the function and the gradient evaluation
  at the same time.  Evaluating both at once may be more efficient that
  evaluating each separately.
*/
PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal* f, Vec G,
                                    void* ptr);

PetscErrorCode FormFunctionGradientTest(Tao tao, Vec X, PetscReal* f, Vec G,
                                        void* ptr);

/*
  FormHessian - Evaluates Hessian matrix.

  Input Parameters:
  .  tao   - the Tao context
  .  x     - input vector
  .  ptr   - optional user-defined context, as set by TaoSetHessian()

  Output Parameters:
  .  H     - Hessian matrix

  Note:  Providing the Hessian may not be necessary.  Only some solvers
  require this matrix.
*/
PetscErrorCode FormHessian(Tao tao, Vec X, Mat H, Mat Hpre, void* ptr);

PetscErrorCode ConvergenceTest(Tao tao, void* ptr);
} // namespace CEPS
