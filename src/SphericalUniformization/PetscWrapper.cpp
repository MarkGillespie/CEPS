#include "PetscWrapper.h"
namespace CEPS {
PetscWrapper::PetscWrapper(PartiallyDecoratedTriangulation& tri_,
                           size_t vInfinityIndex_)
    : tri(tri_), vInfinityIndex(vInfinityIndex_) {
    vInfinity            = tri.mesh->vertex(vInfinityIndex);
    distancesToHorocycle = tri.distanceOfHorocyclesTo(vInfinity);
    setup(std::numeric_limits<double>::infinity());
}

PetscWrapper::PetscWrapper(PartiallyDecoratedTriangulation& tri_,
                           size_t vInfinityIndex_, double vInfinityScaling)
    : tri(tri_), vInfinityIndex(vInfinityIndex_) {
    vInfinity            = tri.mesh->vertex(vInfinityIndex);
    distancesToHorocycle = tri.distanceOfHorocyclesTo(vInfinity);
    setup(vInfinityScaling);
}

PetscWrapper::PetscWrapper(PartiallyDecoratedTriangulation& tri_,
                           VertexData<double>& distancesToHorocycle_,
                           size_t vInfinityIndex_, double vInfinityScaling)
    : tri(tri_), vInfinityIndex(vInfinityIndex_),
      distancesToHorocycle(distancesToHorocycle_) {
    vInfinity = tri.mesh->vertex(vInfinityIndex);
    setup(vInfinityScaling);
}

void PetscWrapper::setup(double vInfinityScaling) {
    vIdx = tri.mesh->getVertexIndices();

    bool needsScaling = false;
    VertexData<double> scaleData(*tri.mesh, 0);
    for (Vertex v : tri.mesh->vertices()) {
        if (v != vInfinity &&
            tri.logScaleFactors[v] <= -distancesToHorocycle[v]) {
            scaleData[v] =
                -distancesToHorocycle[v] - tri.logScaleFactors[v] + 0.2;
            needsScaling = true;
        }
    }

    // if (needsScaling) {
    //     cout << "Mesh was rescaled" << endl;
    //     tri.setScaleFactors(scaleData);
    // } else {
    //     cout << "Mesh was not rescaled" << std::flush << endl;
    // }

    // Scale vInfinity to be infinitely far away
    tri.setVertexScaleFactor(vInfinity, vInfinityScaling);

    if (useFlips) tri.flipToDelaunay(GeometryType::HYPERBOLIC);
    tri.checkEdgeLengthsNaN();
}

bool PetscWrapper::uniformize(int argc, char** argv) {
    // clang-format off
    PetscErrorCode ierr; // used to check for functions returning nonzeros
    PetscReal zero = 0.0;
    Vec x;      // solution vector
    Vec lb, ub; // lower bounds
    Mat H;
    Tao tao; // Tao solver context
    PetscBool flg;
    size_t N = tri.mesh->nVertices();

    static char help[] = "";

    // for (size_t iA =0;(int)iA<argc;++iA) {
    //   std::cout << argv[iA] << std::endl;
    // }

    int size, rank;
    /* Initialize TAO and PETSc */
    PetscInitialize(&argc, &argv, (char*)0, help);

    /* Allocate vectors for the solution, and hessian */
    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x);  CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &lb); CHKERRQ(ierr);

    size_t maxDegree = 0;
    for (Vertex v : tri.mesh->vertices())
      maxDegree = std::max(maxDegree, v.degree());

    // TODO: set number nonzero blocks properly
    ierr = MatCreateSeqBAIJ(PETSC_COMM_SELF, 1, N, N, 2 * maxDegree, NULL, &H); CHKERRQ(ierr);
    ierr = MatSetOption(H, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);


    /* Create TAO solver with desired solution method */
    ierr = TaoCreate(PETSC_COMM_SELF, &tao); CHKERRQ(ierr);
    // ierr = TaoSetType(tao, TAOLMVM);         CHKERRQ(ierr);

    /* Set solution vec and an initial guess */
    ierr = VecSet(x, zero);             CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao, x); CHKERRQ(ierr);

    /* Set routines for function, gradient, hessian evaluation */
    ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient, this); CHKERRQ(ierr);
    ierr = TaoSetHessianRoutine(tao, H, H, FormHessian, this);                 CHKERRQ(ierr);

    /* Set a routine that defines the bounds */
    ierr = VecDuplicate(x, &lb); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &ub); CHKERRQ(ierr);
    PetscReal *lbArr, *ubArr;
    ierr = VecGetArray(lb, &lbArr); CHKERRQ(ierr);
    ierr = VecGetArray(ub, &ubArr); CHKERRQ(ierr);
    for (Vertex v : tri.mesh->vertices()) {
      if (v == vInfinity) {
        lbArr[vIdx[v]] = 0;
        ubArr[vIdx[v]] = 0;
      } else {
        lbArr[vIdx[v]] = -distancesToHorocycle[v];
        ubArr[vIdx[v]] = PETSC_INFINITY;
      }
    }
    ierr = VecRestoreArray(lb, &lbArr);             CHKERRQ(ierr);
    ierr = VecRestoreArray(ub, &ubArr);             CHKERRQ(ierr);
    ierr = TaoSetVariableBounds(tao, lb, ub);       CHKERRQ(ierr);

    // ierr = TaoSetConvergenceTest(tao, ConvergenceTest, this); CHKERRQ(ierr);
    ierr = TaoSetTolerances(tao, 1e-11, 0, 0); CHKERRQ(ierr);
    ierr = TaoSetMaximumIterations(tao, 2000); CHKERRQ(ierr);

    /* Check for TAO command line options */
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

    /* SOLVE THE APPLICATION */
    ierr = TaoSolve(tao); CHKERRQ(ierr);

    /* Check convergence */
    TaoConvergedReason reason;
    ierr = TaoGetConvergedReason(tao, &reason); CHKERRQ(ierr);

    const PetscReal* xArr;

    ierr = VecGetArrayRead(x, &xArr); CHKERRQ(ierr);
    VertexData<double> logScaleFactors(*tri.mesh);
    for (Vertex v : tri.mesh->vertices()) {
      if (v != vInfinity) logScaleFactors[v] = xArr[vIdx[v]];
    }

    logScaleFactors[vInfinity] = tri.logScaleFactors[vInfinity];
    tri.setScaleFactors(logScaleFactors);
    if (useFlips) {
      tri.flipToDelaunay(GeometryType::HYPERBOLIC);
    }
    cleanUpTriangulation();

    ierr = VecRestoreArrayRead(x, &xArr); CHKERRQ(ierr);

    ierr = TaoDestroy(&tao); CHKERRQ(ierr);
    ierr = VecDestroy(&x);   CHKERRQ(ierr);
    ierr = MatDestroy(&H);   CHKERRQ(ierr);

    PetscFinalize();
    return reason == TAO_CONVERGED_GATOL
        || reason == TAO_CONVERGED_GRTOL
        || reason == TAO_CONVERGED_GTTOL
        || reason == TAO_CONVERGED_STEPTOL
        || reason == TAO_CONVERGED_MINF
        || reason == TAO_CONVERGED_USER;
    // clang-format on
}

/*
    FormFunctionGradient - Evaluates the function, f(X), and gradient, G(X).

    Input Parameters:
    .   tao  - the Tao context
    .   X    - input vector
    .   ptr  - optional user-defined context, as set by TaoSetFunctionGradient()

    Output Parameters:
    .   G - vector containing the newly evaluated gradient
    .   f - function value

    Note:
    Some optimization methods ask for the function and the gradient evaluation
    at the same time.  Evaluating both at once may be more efficient that
    evaluating each separately.
*/
PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal* f, Vec G,
                                    void* ptr) {
    // clang-format off
    PetscWrapper* user = (PetscWrapper*)ptr;
    PetscErrorCode ierr;
    const PetscReal* x;
    PetscReal* g;

    /* Get pointers to vector data */
    ierr = VecGetArrayRead(X, &x); CHKERRQ(ierr);
    ierr = VecGetArray(G, &g);     CHKERRQ(ierr);

    VertexData<double> logScaleFactors(*user->tri.mesh);
    for (Vertex v : user->tri.mesh->vertices()) {
      if (v != user->vInfinity) logScaleFactors[v] = x[user->vIdx[v]];
    }
    // for (Vertex v : user->vInfinity.adjacentVertices())
    // {std::cout << "setting " << v << " to " << x[user->vIdx[v]]<<vendl;}

    logScaleFactors[user->vInfinity] =
        user->tri.logScaleFactors[user->vInfinity];
    if (user->useFlips) {
      std::clock_t start = std::clock();
      user->tri.setScaleFactors(logScaleFactors);
      user->tri.flipToDelaunay(GeometryType::HYPERBOLIC);
      double flippingDuration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
      user->flipTime += flippingDuration;
    } else {
      user->tri.setScaleFactors(logScaleFactors);
    }

    /* Compute G(X) */
    Vector<double> gradVec = sphericalGradient(user->tri);
    for (Eigen::Index i = 0; i < gradVec.size(); ++i) {
      g[i] = gradVec[i];
    }
    g[user->vInfinityIndex] = 0;


    /* Restore vectors */
    ierr = VecRestoreArrayRead(X, &x); CHKERRQ(ierr);
    ierr = VecRestoreArray(G, &g);     CHKERRQ(ierr);

    *f = sphericalObjective(user->tri);

    // std::cout << "Objective: " << *f << "\t Gradient norm: " << gradVec.norm() << vendl;

    return 0;
    // clang-format on
}

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

PetscErrorCode FormHessian(Tao tao, Vec X, Mat H, Mat Hpre, void* ptr) {
    // clang-format off
    PetscWrapper* user = (PetscWrapper*)ptr;
    PetscErrorCode ierr;
    const PetscReal* x;
    PetscBool assembled;

    /* Zero existing matrix entries */
    ierr = MatAssembled(H, &assembled); CHKERRQ(ierr);
    if (assembled) {
        ierr = MatZeroEntries(H); CHKERRQ(ierr);
    }

    /* Get a pointer to vector data */
    ierr = VecGetArrayRead(X, &x); CHKERRQ(ierr);

    VertexData<double> logScaleFactors(*user->tri.mesh);
    for (Vertex v : user->tri.mesh->vertices()) {
      if (v != user->vInfinity) logScaleFactors[v] = x[user->vIdx[v]];
    }
    logScaleFactors[user->vInfinity] =
        user->tri.logScaleFactors[user->vInfinity];
    if (user->useFlips) {
      std::clock_t start = std::clock();
      user->tri.setScaleFactors(logScaleFactors);
      user->tri.flipToDelaunay(GeometryType::HYPERBOLIC);
      double flippingDuration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
      user->flipTime += flippingDuration;
    } else {
      user->tri.setScaleFactors(logScaleFactors);
    }

    user->adjustTriangulation();
    SparseMatrix<double> h = sphericalHessian(user->tri);

    /* Compute H(X) entries */
    for (Eigen::Index k = 0; k < h.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(h, k); it; ++it) {
            size_t row   = it.row();
            size_t col   = it.col();
            double value = it.value();
            ierr         = MatSetValue(H, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    int iV = user->vInfinityIndex;
    ierr = MatSetValue(H, iV, iV, 1, INSERT_VALUES); CHKERRQ(ierr);


    ierr = VecRestoreArrayRead(X, &x); CHKERRQ(ierr);

    /* Assemble matrix */
    ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);

    return 0;
    // clang-format on
}

PetscErrorCode PetscWrapper::testOptimize(int argc, char** argv) {
    PetscErrorCode ierr; // used to check for functions returning nonzeros
    PetscReal zero = 0.0;
    Vec x; // solution vector
    PetscReal* xArr;
    Vec lb, ub; // lower bounds
    Mat H;
    Tao tao; // Tao solver context
    PetscBool flg;
    size_t N = 3;

    static char help[] =
        "This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on a single processor.  We \n\
minimize the extended Rosenbrock function: \n\
   sum_{i=0}^{n/2-1} ( alpha*(x_{2i+1}-x_{2i}^2)^2 + (1-x_{2i})^2 ) \n";


    int size, rank;
    /* Initialize TAO and PETSc */
    PetscInitialize(&argc, &argv, (char*)0, help);

    /* Allocate vectors for the solution, and hessian */
    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x);
    CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, N, &lb);
    CHKERRQ(ierr);
    ierr = MatCreateSeqBAIJ(PETSC_COMM_SELF, 1, N, N, 2, NULL, &H);
    CHKERRQ(ierr);

    /* The TAO code begins here */

    /* Create TAO solver with desired solution method */
    ierr = TaoCreate(PETSC_COMM_SELF, &tao);
    CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOLMVM);
    CHKERRQ(ierr);

    /* Set solution vec and an initial guess */
    ierr = VecGetArray(x, &xArr);
    for (size_t i = 0; i < N; ++i) {
        xArr[i] = 8;
    }
    ierr = VecRestoreArray(x, &xArr);
    CHKERRQ(ierr);

    ierr = TaoSetInitialVector(tao, x);
    CHKERRQ(ierr);

    /* Set routines for function, gradient, hessian evaluation */
    ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradientTest, &N);
    CHKERRQ(ierr);

    /* Set a routine that defines the bounds */
    ierr = VecDuplicate(x, &lb);
    ierr = VecDuplicate(x, &ub);
    CHKERRQ(ierr);
    PetscReal *lbArr, *ubArr;
    ierr = VecGetArray(lb, &lbArr);
    ierr = VecGetArray(ub, &ubArr);
    CHKERRQ(ierr);
    for (size_t i = 0; i < N; ++i) {
        lbArr[i] = 2;
        ubArr[i] = PETSC_INFINITY;
    }
    ierr = VecRestoreArray(lb, &lbArr);
    ierr = VecRestoreArray(ub, &ubArr);
    CHKERRQ(ierr);
    ierr = TaoSetVariableBounds(tao, lb, ub);
    CHKERRQ(ierr);

    ierr = TaoSetTolerances(tao, 1e-12, 1e-12, 1e-12);
    CHKERRQ(ierr);

    /* Check for TAO command line options */
    ierr = TaoSetFromOptions(tao);
    CHKERRQ(ierr);

    /* SOLVE THE APPLICATION */
    ierr = TaoSolve(tao);
    CHKERRQ(ierr);

    ierr = VecGetArray(x, &xArr);
    for (size_t i = 0; i < N; ++i) {
        cout << xArr[i] << "\t";
    }
    cout << endl;
    ierr = VecRestoreArray(x, &xArr);
    CHKERRQ(ierr);


    ierr = TaoDestroy(&tao);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = MatDestroy(&H);
    CHKERRQ(ierr);

    PetscFinalize();
    cout << "PETSc DONE" << endl;
    return 0;
}

PetscErrorCode FormFunctionGradientTest(Tao tao, Vec X, PetscReal* f, Vec G,
                                        void* ptr) {
    int N = *((int*)ptr);
    PetscErrorCode ierr;
    const PetscReal* x;
    PetscReal* g;

    /* Get pointers to vector data */
    ierr = VecGetArrayRead(X, &x);
    CHKERRQ(ierr);
    ierr = VecGetArray(G, &g);
    CHKERRQ(ierr);

    /* Compute G(X) */
    *f = 0;
    for (size_t i = 0; (int)i < N; ++i) {
        *f += x[i] * x[i];
        g[i] = 2 * x[i];
        cout << x[i] << "\t";
    }
    cout << N << endl;

    /* Restore vectors */
    ierr = VecRestoreArrayRead(X, &x);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(G, &g);
    CHKERRQ(ierr);

    return 0;
}

PetscErrorCode ConvergenceTest(Tao tao, void* ptr) {
    // clang-format off
    PetscErrorCode ierr;
    PetscReal gatol, grtol, gttol, f, gnorm, cnorm, xdiff;
    PetscInt its, maxits;
    TaoConvergedReason reason;
    PetscWrapper* user = (PetscWrapper*)ptr;

    ierr = TaoGetSolutionStatus(tao, &its, &f, &gnorm, &cnorm, &xdiff, &reason);
    ierr = TaoGetTolerances(tao, &gatol, &grtol, &gttol); CHKERRQ(ierr);
    ierr = TaoGetMaximumIterations(tao, &maxits);

    if (its >= maxits) {
      ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS); CHKERRQ(ierr);
    } else {
      double di, cs, angle;
      user->kktError(di, cs, angle);

      if ((di < gatol && cs < gatol) || its >= maxits)  {
        ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_USER); CHKERRQ(ierr);
      }
    }

    return 0;
    // clang-format on
}

void PetscWrapper::kktError(double& dualInfeasibility,
                            double& complementarySlacknessError,
                            double& angleErr) {
    throw_verbose_runtime_error("Not implemented yet");
    // VertexData<double> gradF = tri.objectiveGradientVertexData();

    // dualInfeasibility           = 0;
    // complementarySlacknessError = 0;
    // angleErr                    = 0;
    // VertexData<double> angleSum = tri.angleSums();
    // for (Vertex v : tri.mesh->vertices()) {
    //     if (v != vInfinity) {
    //         double delta = distancesToHorocycle[v];
    //         double u     = tri.logScaleFactors[v];
    //         double g     = -delta - u;
    //         // Note that grad g is a Kronecker delta function So we know that
    //         mu
    //         // is just the corresponding component of grad f
    //         double mu = gradF[v];
    //         // my_assert(g <= 0, "point is infeasible");

    //         if (!adjacent(v, vInfinity)) {
    //             double explicitAngleViolation = 2 * PI - angleSum[v];
    //             // my_assert(
    //             //     abs(explicitAngleViolation - gradF[v]) < 1e-10,
    //             //     "angle violation wrong? diff is " +
    //             //         std::to_string(abs(explicitAngleViolation -
    //             //         gradF[v])));
    //             angleErr = fmax(angleErr, explicitAngleViolation);
    //         }

    //         if (mu < 0) dualInfeasibility = fmax(dualInfeasibility, abs(mu));
    //         complementarySlacknessError =
    //             fmax(abs(mu * g), complementarySlacknessError);
    //     }
    // }
}

bool PetscWrapper::activeConstraint(Vertex v) {
    return abs(tri.logScaleFactors[v] + distancesToHorocycle[v]) < 1e-8;
}

void PetscWrapper::adjustTriangulation() {
    bool done;
    do {
        done = true;
        for (Halfedge spoke : vInfinity.outgoingHalfedges()) {
            Halfedge he = spoke.next();
            if (activeConstraint(he.vertex()) &&
                activeConstraint(he.next().vertex()) &&
                !tri.isEssential(he.edge())) {
                tri.flipEdge(he.edge(), GeometryType::HYPERBOLIC);
                done = false;
                break;
            }
        }

    } while (!done);
}

// Run after running Petsc optimizer. Petsc doesn't always give me a solution
// that's exactly on the boundary. So I pin the boundary vertices to our
// constraints. This can mess up the angle sums a little bit, so I run a few
// rounds of ordinary CETM with fixed boundary conditions to clean everything up
void PetscWrapper::cleanUpTriangulation() {
    double di, cs, angle;
    adjustTriangulation();

    double maxClamp = 0;
    // This seems to mess up angles?
    // Clamp boundary vertices to constraints
    for (Halfedge he : vInfinity.incomingHalfedges()) {
        Vertex v = he.vertex();
        maxClamp = fmax(maxClamp,
                        abs(tri.logScaleFactors[v] + distancesToHorocycle[v]));
        tri.setVertexScaleFactor(v, -distancesToHorocycle[v]);
    }
    // cout << "maxClamp: " << maxClamp << endl;
    // tri.setVertexScaleFactor(vInfinity, 1);
    tri.updateEdgeLengths();

    // kktError(di, cs, angle);
    // if (angle > 1e-8) {
    //   cerr << "angle err: " << angle << endl;
    // }

    // my_assert(tri.satisfiesTriangleInequality(),
    //           "violates triangle inequality");
}
} // namespace CEPS
