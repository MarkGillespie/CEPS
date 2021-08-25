#pragma once

#include <gtest/gtest.h>

#include <Eigen/SparseCore>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <stdlib.h> /* rand */

#define EXPECT_VEC2_EQ(a, b) EXPECT_PRED2(Vector2Eq, a, b)
#define EXPECT_VEC3_EQ(a, b) EXPECT_PRED2(Vector3Eq, a, b)
#define EXPECT_VEC2_NEAR(a, b, tol) EXPECT_PRED3(Vector2Near, a, b, tol)
#define EXPECT_VEC3_NEAR(a, b, tol) EXPECT_PRED3(Vector3Near, a, b, tol)
#define EXPECT_MAT_EQ(a, b) ExpectMatEq(a, b)
#define EXPECT_MAT_NEAR(a, b, tol) ExpectMatEq(a, b, tol)
#define ASSERT_MAT_FINITE(m) AssertMatFinite(m)

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

double fRand(double fMin, double fMax) {
    double f = (double)rand() / (RAND_MAX + 1.0);
    return fMin + f * (fMax - fMin);
}

// =============================================================
//                Template Magic for Eigen
// =============================================================
// googletest doesn't like printing out eigen matrices when they fail tests
// so I stole this code from
// https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix

template <class Base>
class EigenPrintWrap : public Base {
    friend void PrintTo(const EigenPrintWrap& m, ::std::ostream* o) {
        size_t width  = (m.cols() < 10) ? m.cols() : 10;
        size_t height = (m.rows() < 10) ? m.rows() : 10;
        *o << "\n" << m.topLeftCorner(height, width) << "...";
    }
};

template <class Base>
const EigenPrintWrap<Base>& print_wrap(const Base& base) {
    return static_cast<const EigenPrintWrap<Base>&>(base);
}

// TODO: upgrade to work on arbitrary matrix types (or at least integer
// matrices)
bool MatrixEq(const EigenPrintWrap<Eigen::MatrixXd>& lhs_,
              const EigenPrintWrap<Eigen::MatrixXd>& rhs_, double difference,
              double threshold = -1) {
    Eigen::MatrixXd lhs = static_cast<Eigen::MatrixXd>(lhs_);
    Eigen::MatrixXd rhs = static_cast<Eigen::MatrixXd>(rhs_);
    double err          = (lhs - rhs).norm();
    if (threshold > 0) {
        bool equal = abs(err) < threshold;
        if (!equal) cerr << "norm of difference: " << err << endl;
        return equal;
    } else {
        const ::testing::internal::FloatingPoint<double> difference(err),
            zero(0);
        bool equal = difference.AlmostEquals(zero);
        if (!equal) cerr << "norm of difference: " << err << endl;
        return equal;
    }
}

void ExpectMatEq(const Eigen::MatrixXd& a_, const Eigen::MatrixXd& b_,
                 double threshold = -1) {
    EXPECT_PRED4(MatrixEq, print_wrap(a_), print_wrap(b_), (a_ - b_).norm(),
                 threshold);
}

// Checks that matrix is finite and not nan
bool MatFinite(const Eigen::MatrixXd& m) {
    return std::isfinite(m.squaredNorm());
}

void AssertMatFinite(const Eigen::MatrixXd& m_) {
    ASSERT_PRED1(MatFinite, print_wrap(m_));
}

// =============================================================
//    Floating Point Comparison for Geometry-Central Vectors
// =============================================================

// Vector2 floating point equality (up to sign)
bool Vector2EqPm(const Vector2& a, const Vector2& b) {
    const ::testing::internal::FloatingPoint<double> ax(a.x), ay(a.y), bx(b.x),
        by(b.y), nbx(-b.x), nby(-b.y);
    return (ax.AlmostEquals(bx) && ay.AlmostEquals(by)) ||
           (ax.AlmostEquals(nbx) && ay.AlmostEquals(nby));
}

// Vector2 floating point equality
bool Vector2Eq(const Vector2& a, const Vector2& b) {
    const ::testing::internal::FloatingPoint<double> ax(a.x), ay(a.y), bx(b.x),
        by(b.y);
    return ax.AlmostEquals(bx) && ay.AlmostEquals(by);
}

// Vector2 floating point near
bool Vector2Near(const Vector2& a, const Vector2& b, const double& tol) {
    return (a - b).norm() < tol;
}

// Vector3 floating point equality
bool Vector3Eq(const Vector3& a, const Vector3& b) {
    const ::testing::internal::FloatingPoint<double> ax(a.x), ay(a.y), bx(b.x),
        by(b.y), az(a.z), bz(b.z);
    return ax.AlmostEquals(bx) && ay.AlmostEquals(by) && az.AlmostEquals(bz);
}

// Vector3 floating point near
bool Vector3Near(const Vector3& a, const Vector3& b, const double& tol) {
    return (a - b).norm() < tol;
}
