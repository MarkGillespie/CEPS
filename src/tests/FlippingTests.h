#include "test_utils.h"

#include "NumericalNormalCoordinatesComputations.h"
#include "Tracing.h"
#include "Triangulation.h"
#include "Utils.h"

#include "geometrycentral/surface/meshio.h"

#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"

using namespace CEPS;

class FlippingTest : public ::testing::Test {
  public:
    static std::unique_ptr<Triangulation> tri;
    static std::unique_ptr<VertexData<Vector3>> vertexPositions;

  protected:
    static void SetUpTestSuite() {
        auto example    = loadMesh("../src/tests/bunny_tiny.obj");
        tri             = std::make_unique<Triangulation>(*std::get<0>(example),
                                              *std::get<1>(example));
        vertexPositions = std::make_unique<VertexData<Vector3>>(
            *tri->mesh, std::get<1>(example)->inputVertexPositions.raw());
    }
};

std::unique_ptr<Triangulation> FlippingTest::tri                   = nullptr;
std::unique_ptr<VertexData<Vector3>> FlippingTest::vertexPositions = nullptr;

TEST_F(FlippingTest, FlippingTwiceJustReversesEdge) {
    Edge e                                = tri->mesh->edge(12);
    EdgeData<double> oldEdgeLengths       = tri->originalEdgeLengths;
    EdgeData<size_t> oldNormalCoordinates = tri->normalCoordinates;
    HalfedgeData<size_t> oldRoundabouts   = tri->roundaboutIndices;

    for (int i = 0; i < 2; ++i) {
        tri->flipEdge(e, GeometryType::EUCLIDEAN);
    }

    // Reversing the edge orientation switched the roundabouts
    std::swap(oldRoundabouts[e.halfedge()],
              oldRoundabouts[e.halfedge().twin()]);

    EXPECT_MAT_NEAR(oldEdgeLengths.toVector().cast<double>(),
                    tri->originalEdgeLengths.toVector().cast<double>(), 1e-8);
    EXPECT_MAT_EQ(oldNormalCoordinates.toVector().cast<double>(),
                  tri->normalCoordinates.toVector().cast<double>());
    EXPECT_MAT_EQ(oldRoundabouts.toVector().cast<double>(),
                  tri->roundaboutIndices.toVector().cast<double>());
}

TEST_F(FlippingTest, EuclideanFlipEdgeLength) {
    Vector2 qi{0, 0};
    Vector2 qj{1, 0};
    Vector2 qk{0.5, 1};
    Vector2 ql{0.5, -1};

    double lij = (qi - qj).norm();
    double ljk = (qk - qj).norm();
    double lki = (qi - qk).norm();
    double lil = (ql - qi).norm();
    double llj = (qj - ql).norm();
    std::array<double, 5> lengths{lij, ljk, lki, lil, llj};

    double flippedLen = flipEuclideanLength(lengths);

    EXPECT_DOUBLE_EQ(flippedLen, 2);
}

TEST_F(FlippingTest, RandomEuclideanFlipEdgeLength) {
    for (size_t i = 0; i < 10; ++i) {
        Vector2 qi{0, 0};
        Vector2 qj{fRand(-1, 1), fRand(-1, 1)};
        Vector2 qk{fRand(-1, 1), fRand(-1, 1)};
        Vector2 ql{fRand(-1, 1), fRand(-1, 1)};

        // Make sure ql is on the opposite side of the line qi-qj as qk is
        double c1 = cross(qj, qk);
        double c2 = cross(qj, ql);
        if (c1 * c2 > 0) {
            // Reflect ql across line qj
            // n is the normal vector to the line;
            Vector2 n{-qj.y, qj.x};
            n = n.normalize();

            ql = ql - 2 * dot(n, ql) * n;
        }

        double lij = (qi - qj).norm();
        double ljk = (qk - qj).norm();
        double lki = (qi - qk).norm();
        double lil = (ql - qi).norm();
        double llj = (qj - ql).norm();
        std::array<double, 5> lengths{lij, ljk, lki, lil, llj};

        double flippedLen = flipEuclideanLength(lengths);

        if (flippedLen > 0) {
            EXPECT_NEAR(flippedLen, (qk - ql).norm(), 1e-8);
        }
    }
}


TEST_F(FlippingTest, SimpleNormalCoordinateExample1) {
    size_t nij = 4;
    size_t njk = 3;
    size_t nki = 5;
    size_t nil = 1;
    size_t nlj = 0;

    size_t nlk = flipNormalCoordinate({nij, njk, nki, nil, nlj});

    ASSERT_EQ(nlk, 2);
}

TEST_F(FlippingTest, SimpleNormalCoordinateExample2) {
    // Tricky because 1 edge goes from l to k
    size_t nij = 49;
    size_t njk = 1;
    size_t nki = 40;
    size_t nil = 19;
    size_t nlj = 8;

    size_t nlk = flipNormalCoordinate({nij, njk, nki, nil, nlj});

    ASSERT_EQ(nlk, 0);
}

TEST_F(FlippingTest, RandomNormalCoordinates) {
    srand(0); // use the same seed every time the test is run

    size_t maxN = 50;
    for (size_t i = 0; i < 1000; ++i) {
        size_t a = rand() % maxN;
        size_t b = rand() % maxN;
        size_t c = rand() % maxN;
        size_t d = rand() % maxN;
        size_t e = rand() % maxN;

        // If abe satisfies the triangle inequality, then there should be an
        // even number of crossings (i.e. a + b + e is even) since every
        // line crosses two edges. If our random numbers satisfy the
        // triangle inequality but sum to an odd number, we can just add one
        // to the minimum. Note that if we satisfy the triangle inequality
        // with equality, the sum is even, so we don't have to worry about
        // accidentally causing a triangle inequality violation
        if (triangleIneq(a, b, e) && (a + b + e) % 2 == 1) {
            if (a <= b) {
                ++a;
            } else {
                ++b;
            }
        }

        if (triangleIneq(c, d, e) && (c + d + e) % 2 == 1) {
            if (c <= d) {
                ++c;
            } else {
                ++d;
            }
        }

        testCoords(a, b, c, d, e);
    }
}

TEST_F(FlippingTest, NormalCoordinatesAgreeWithRoundaboutDegrees) {
    tri->flipToDelaunay(GeometryType::EUCLIDEAN);
    auto em = [&](Corner c) -> size_t {
        return positivePart(
            tri->normalCoordinates[c.halfedge().next().edge()] -
            tri->normalCoordinates[c.halfedge().next().next().edge()] -
            tri->normalCoordinates[c.halfedge().edge()]);
    };
    for (Vertex v : tri->mesh->vertices()) {
        size_t normalCoordinateDegree = 0;
        for (Halfedge he : v.outgoingHalfedges()) {
            normalCoordinateDegree += em(he.corner());
            if (tri->normalCoordinates[he.edge()] == 0)
                normalCoordinateDegree++;
        }
        EXPECT_EQ(normalCoordinateDegree, tri->roundaboutDegrees[v]);
    }
}
