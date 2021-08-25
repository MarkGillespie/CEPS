#include "test_utils.h"

#include "CommonRefinement.h"

using namespace CEPS;
using namespace CEPS::ImplementationDetails;

class CommonRefinementTests : public testing::Test {};

TEST_F(CommonRefinementTests, doubleRefineTemplatingIsOkay) {

    std::array<std::vector<double>, 3> aLines{
        std::vector<double>{0, 0.2, 0.4, 1},
        std::vector<double>{0, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1},
        std::vector<double>{0, 0.1, 0.8, 1}};
    std::array<std::vector<double>, 3> bLines{
        std::vector<double>{0, 0.3, 0.6, 0.9, 1},
        std::vector<double>{0, 0.1, 0.7, 1},
        std::vector<double>{0, 0.2, 0.4, 0.5, 1},
    };

    std::array<std::vector<std::array<Vector2, 2>>, 3> val1;
    std::array<std::vector<std::array<Vector3, 2>>, 3> val2;

    for (size_t iE = 0; iE < 3; ++iE) {
        val1[iE] = std::vector<std::array<Vector2, 2>>(aLines[iE].size());
        val2[iE] = std::vector<std::array<Vector3, 2>>(bLines[iE].size());
    }

    Vector2 vZero{0, 0};
    std::vector<Vector2> longSideCornerVal1{vZero, vZero, vZero};
    std::vector<Vector3> longSideCornerVal2{};


    size_t nVertices;
    std::vector<std::vector<size_t>> polygons;
    std::vector<std::vector<Vector2>> cornerVal1;
    std::vector<std::vector<Vector3>> cornerVal2;
    std::tie(nVertices, polygons, cornerVal1, cornerVal2) = sliceTri(
        aLines, bLines, val1, val2, longSideCornerVal1, longSideCornerVal2);
}

TEST_F(CommonRefinementTests, splitTriangle) {

    std::array<std::vector<double>, 3> aLines{
        std::vector<double>{0, 0.2, 0.4, 1},
        std::vector<double>{0, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1},
        std::vector<double>{0, 0.1, 0.8, 1}};
    std::array<std::vector<double>, 3> bLines{
        std::vector<double>{0, 0.3, 0.6, 0.9, 1},
        std::vector<double>{0, 0.1, 0.7, 1},
        std::vector<double>{0, 0.2, 0.4, 0.5, 1},
    };

    std::array<std::vector<std::array<double, 2>>, 3> val1, val2;
    for (size_t i = 0; i < 3; ++i) {
        val1[i].reserve(aLines[i].size());
        for (auto l : aLines[i]) {
            val1[i].push_back(std::array<double, 2>{0, 0});
        }
        val2[i].reserve(bLines[i].size());
        for (auto l : bLines[i]) {
            val2[i].push_back(std::array<double, 2>{0, 0});
        }
    }

    std::vector<double> longSideCornerVal1{0, 0, 0};
    std::vector<double> longSideCornerVal2{};


    size_t nVertices;
    std::vector<std::vector<size_t>> polygons;
    std::vector<std::vector<double>> cornerVal1, cornerVal2;
    std::tie(nVertices, polygons, cornerVal1, cornerVal2) = sliceTri(
        aLines, bLines, val1, val2, longSideCornerVal1, longSideCornerVal2);

    // Computed with a working version of the function
    ASSERT_EQ(nVertices, 31);
    ASSERT_EQ(polygons.size(), 21);

    auto distinct = [](const std::vector<size_t>& a,
                       const std::vector<size_t>& b) {
        if (a.size() != b.size()) {
            return true;
        } else {
            for (size_t i : a) {
                bool bContainsI = false;
                for (size_t j : b) {
                    if (i == j) bContainsI = true;
                }
                if (!bContainsI) {
                    return true;
                }
            }
            return false;
        }
    };

    // Check that polygons are unique
    for (size_t iP = 0; iP < polygons.size(); ++iP) {
        std::vector<size_t> p = polygons[iP];
        for (size_t iQ = 0; iQ < iP; ++iQ) {
            std::vector<size_t> q = polygons[iQ];
            ASSERT_TRUE(distinct(p, q));
        }
    }
}
