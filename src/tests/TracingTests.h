#include "test_utils.h"

#include "Tracing.h"
#include "Triangulation.h"
#include "Utils.h"

#include "geometrycentral/surface/meshio.h"

#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"

using namespace CEPS;

class TracingTest : public ::testing::Test {
  public:
    static std::unique_ptr<Triangulation> T1;
    static std::unique_ptr<Triangulation> T2;
    static std::unique_ptr<VertexData<Vector3>> T1Positions;
    double localDouble;

  protected:
    static void SetUpTestSuite() {
        auto example = loadMesh("../src/tests/bunny_tiny.obj");
        T1           = std::make_unique<Triangulation>(*std::get<0>(example),
                                             *std::get<1>(example));
        T2           = std::make_unique<Triangulation>(*T1);
        T2->flipToDelaunay(GeometryType::EUCLIDEAN);
        T1Positions = std::make_unique<VertexData<Vector3>>(
            *T2->mesh, std::get<1>(example)->inputVertexPositions.raw());
    }

    void SetUp() override { localDouble = 10; }
};

std::unique_ptr<Triangulation> TracingTest::T1                = nullptr;
std::unique_ptr<Triangulation> TracingTest::T2                = nullptr;
std::unique_ptr<VertexData<Vector3>> TracingTest::T1Positions = nullptr;

TEST_F(TracingTest, TrivialPathsGetTraced) {
    EdgeData<MeshPath> paths = traceTopologicalTriangulation(*T1, *T1);
    for (Edge e : T1->mesh->edges()) {
        ASSERT_TRUE(paths[e].start.corner == e.halfedge().corner());
        ASSERT_TRUE(paths[e].points.empty());
    }
}

TEST_F(TracingTest, PathsGetTraced) {
    EdgeData<MeshPath> paths = traceTopologicalTriangulation(*T1, *T2);
    for (Edge e : T2->mesh->edges()) {
        ASSERT_TRUE(paths[e].start.corner != Corner());
    }
}

TEST_F(TracingTest, PathsIntersectHalfedgesPositively) {

    EdgeData<MeshPath> paths = traceTopologicalTriangulation(*T1, *T2);
    for (Edge e : T2->mesh->edges()) {
        for (size_t iH = 0; iH + 1 < paths[e].points.size(); ++iH) {
            Halfedge heCurr = paths[e].points[iH].halfedge;
            Halfedge heNext = paths[e].points[iH + 1].halfedge;
            EXPECT_TRUE(heNext == heCurr.twin().next() ||
                        heNext == heCurr.twin().next().next());
        }
    }
}

TEST_F(TracingTest, TransposedPathsIntersectHalfedgesPositively) {

    EdgeData<MeshPath> paths = traceTopologicalTriangulation(*T1, *T2);
    EdgeData<MeshPath> transposedPaths =
        ImplementationDetails::transpose(paths, *T1, *T2);
    for (Edge e : T1->mesh->edges()) {
        for (size_t iH = 0; iH + 1 < transposedPaths[e].points.size(); ++iH) {
            Halfedge heCurr = transposedPaths[e].points[iH].halfedge;
            Halfedge heNext = transposedPaths[e].points[iH + 1].halfedge;
            EXPECT_TRUE(heNext == heCurr.twin().next() ||
                        heNext == heCurr.twin().next().next());
        }
    }
}

TEST_F(TracingTest, TracedEdgesHaveCorrectNumberOfIntersections) {

    EdgeData<MeshPath> T1overT2 = traceTopologicalTriangulation(*T1, *T2);
    EdgeData<MeshPath> T2overT1 =
        ImplementationDetails::transpose(T1overT2, *T1, *T2);

    EdgeData<MeshPath> paths = traceTopologicalTriangulation(*T1, *T2);
    EdgeData<size_t> tracedIntersections(*T2->mesh, 0);
    for (Edge e : T1->mesh->edges()) {
        for (HalfedgePt hPt : paths[e].points) {
            tracedIntersections[hPt.halfedge.edge()]++;
        }
    }
    for (Edge e : T2->mesh->edges()) {
        EXPECT_EQ(tracedIntersections[e], T2->normalCoordinates[e]);
    }
}
