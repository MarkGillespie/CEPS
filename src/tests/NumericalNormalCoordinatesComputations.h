#pragma once
#include "FlipFormulas.h"
#include "Utils.h"

using namespace CEPS;

// ============================================
//        Coordinate Geometry Tests
// ============================================

bool triangleIneq(size_t a, size_t b, size_t e) {
    return (a + b >= e && b + e >= a && e + a >= b);
}

// a, b, e should be in counterclockwise order
void divideEdges(size_t a, size_t b, size_t e, size_t& a2b, size_t& b2e,
                 size_t& e2a, size_t& amid, size_t& bmid, size_t& emid) {
    if (triangleIneq(a, b, e)) {
        if ((a + b + e) % 2 == 1) {
            cout << "a: " << a << "\tb: " << b << "\te: " << e << endl;
        }
        assert((a + b - e) % 2 == 0);
        assert((b + e - a) % 2 == 0);
        assert((e + a - b) % 2 == 0);

        a2b  = (a + b - e) / 2;
        b2e  = (b + e - a) / 2;
        e2a  = (e + a - b) / 2;
        amid = 0;
        bmid = 0;
        emid = 0;
    } else if (e > a + b) {
        // If extra edges go through e, there can be no edges from a to b
        a2b  = 0;
        b2e  = b;
        e2a  = a;
        amid = 0;
        bmid = 0;
        emid = e - a - b;
    } else if (a > b + e) {
        a2b  = b;
        b2e  = 0;
        e2a  = e;
        amid = a - b - e;
        bmid = 0;
        emid = 0;
    } else { // b > a + e
        assert(b > a + e);
        a2b  = a;
        b2e  = e;
        e2a  = 0;
        amid = 0;
        bmid = b - a - e;
        emid = 0;
    }
    assert(b2e + emid + e2a == e);
}

Vector2 defaultX{1, 0};
Vector2 defaultY{0, 1};
Vector2 defaultZ = {-1, 0};
Vector2 defaultW = {0, -1};
size_t computeNumericalEdgeIntersectionNumber(
    size_t a, size_t b, size_t c, size_t d, size_t e, Vector2 x = defaultX,
    Vector2 y = defaultY, Vector2 z = defaultZ, Vector2 w = defaultW) {

    /*
     * We lay out the triangles so that edge f lies along the y axis. Then
     * we can count intersections with f by counting how many segments cross
     * the y-axis. Note that this means we can do the computation separately
     * in the two triangles.
     *
     *                 y =(0, 1)
     *                 /\
     *                /  \
     *               /    \
     *              b      a
     *             /   U    \
     *            /          \
     * (-1, 0) = z ---- e ---- x = (1, 0)
     *            \          /
     *             \   L    /
     *              c      d
     *               \    /
     *                \  /
     *                 \/
     *                 w = (0, -1)
     */
    // Indices assume clockwise orientation
    char aSegmentEndpointEdge[a];
    size_t aSegmentEndpointIndex[a];
    char bSegmentEndpointEdge[b];
    size_t bSegmentEndpointIndex[b];

    char cSegmentEndpointEdge[c];
    size_t cSegmentEndpointIndex[c];
    char dSegmentEndpointEdge[d];
    size_t dSegmentEndpointIndex[d];

    char eSegmentEndpointEdge[e]; // segments in lower triangle
    size_t eSegmentEndpointIndex[e];

    size_t a2b, b2e, e2a, amid, bmid, emid;
    divideEdges(a, b, e, a2b, b2e, e2a, amid, bmid, emid);

    for (size_t i = 0; i < a2b; ++i) {
        aSegmentEndpointEdge[e2a + amid + i]  = 'b';
        aSegmentEndpointIndex[e2a + amid + i] = i;

        bSegmentEndpointEdge[i]  = 'a';
        bSegmentEndpointIndex[i] = e2a + amid + i;
    }

    for (size_t i = 0; i < b2e; ++i) {
        bSegmentEndpointEdge[a2b + bmid + i]  = 'e';
        bSegmentEndpointIndex[a2b + bmid + i] = i;
    }

    for (size_t i = 0; i < e2a; ++i) {
        aSegmentEndpointEdge[i]  = 'e';
        aSegmentEndpointIndex[i] = b2e + emid + i;
        assert(b2e + emid + i < e);
    }

    for (size_t i = 0; i < amid; ++i) {
        aSegmentEndpointEdge[e2a + i]  = 'z';
        aSegmentEndpointIndex[e2a + i] = 0;
    }

    for (size_t i = 0; i < bmid; ++i) {
        bSegmentEndpointEdge[a2b + i]  = 'x';
        bSegmentEndpointIndex[a2b + i] = 0;
    }

    size_t c2d, d2e, e2c, cmid, dmid, emidDown;
    divideEdges(c, d, e, c2d, d2e, e2c, cmid, dmid, emidDown);

    for (size_t i = 0; i < c2d; ++i) {
        cSegmentEndpointEdge[e2c + cmid + i]  = 'd';
        cSegmentEndpointIndex[e2c + cmid + i] = i;

        dSegmentEndpointEdge[i]  = 'c';
        dSegmentEndpointIndex[i] = e2c + cmid + i;
    }

    for (size_t i = 0; i < d2e; ++i) {
        dSegmentEndpointEdge[c2d + dmid + i]  = 'e';
        dSegmentEndpointIndex[c2d + dmid + i] = i;

        eSegmentEndpointEdge[e - 1 - i]  = 'd';
        eSegmentEndpointIndex[e - 1 - i] = c2d + dmid + i;
    }

    for (size_t i = 0; i < e2c; ++i) {
        cSegmentEndpointEdge[i]  = 'e';
        cSegmentEndpointIndex[i] = d2e + emidDown + i;

        eSegmentEndpointEdge[i]  = 'c';
        eSegmentEndpointIndex[i] = i;
    }

    for (size_t i = 0; i < cmid; ++i) {
        cSegmentEndpointEdge[e2c + i]  = 'x';
        cSegmentEndpointIndex[e2c + i] = 0;
    }
    for (size_t i = 0; i < dmid; ++i) {
        dSegmentEndpointEdge[c2d + i]  = 'z';
        dSegmentEndpointIndex[c2d + i] = 0;
    }
    for (size_t i = 0; i < emidDown; ++i) {
        eSegmentEndpointEdge[e2c + i]  = 'w';
        eSegmentEndpointIndex[e2c + i] = 0;
    }

    for (size_t i = 0; i < e; ++i) {
        char end = eSegmentEndpointEdge[i];
        assert(end == 'c' || end == 'w' || end == 'd');
    }

    assert(e2c + emidDown + d2e == e);

    std::vector<std::pair<Vector2, Vector2>> lines;
    auto lineEndpoint = [&](char endpointEdge, size_t endpointIndex) {
        if (endpointEdge == 'x') {
            return x;
        } else if (endpointEdge == 'y') {
            return y;
        } else if (endpointEdge == 'z') {
            return z;
        } else if (endpointEdge == 'w') {
            return w;
        } else if (endpointEdge == 'b') {
            return y + (double)(endpointIndex + 1) / ((double)b + 1) * (z - y);
        } else if (endpointEdge == 'c') {
            return z + (double)(endpointIndex + 1) / ((double)c + 1) * (w - z);
        } else if (endpointEdge == 'd') {
            return w + (double)(endpointIndex + 1) / ((double)d + 1) * (x - w);
        } else if (endpointEdge == 'e') {
            assert(endpointIndex < e);
            if (eSegmentEndpointEdge[endpointIndex] == 'c') {
                return z + (double)(eSegmentEndpointIndex[endpointIndex] + 1) /
                               ((double)c + 1) * (w - z);
            } else if (eSegmentEndpointEdge[endpointIndex] == 'w') {
                return w;
            } else if (eSegmentEndpointEdge[endpointIndex] == 'd') {
                return w + (double)(eSegmentEndpointIndex[endpointIndex] + 1) /
                               ((double)d + 1) * (x - w);
            } else {
                // lines from e downwards should only go to c, the opposite
                // vertex, or d
                cout << "Error: edge from e goes to '"
                     << eSegmentEndpointEdge[endpointIndex] << "'" << endl;
                exit(1);
            }
        } else {
            exit(1);
        }
    };

    // add edges adjacent to a;
    for (size_t i = 0; i < a; ++i) {
        Vector2 start = x + (double)(i + 1) / ((double)a + 1) * (y - x);
        Vector2 end =
            lineEndpoint(aSegmentEndpointEdge[i], aSegmentEndpointIndex[i]);
        lines.emplace_back(std::make_pair(start, end));
    }
    // add edges adjacent to b (which don't go to a)
    for (size_t i = 0; i < b; ++i) {
        Vector2 start = y + (double)(i + 1) / ((double)b + 1) * (z - y);
        if (bSegmentEndpointEdge[i] == 'a') {
            continue;
        }
        Vector2 end =
            lineEndpoint(bSegmentEndpointEdge[i], bSegmentEndpointIndex[i]);
        lines.emplace_back(std::make_pair(start, end));
    }

    // add edges coming out of vertex y (these pass through e at positions
    // b2e to b2e + emid)
    for (size_t j = 0; j < emid; ++j) {
        Vector2 start = y;
        Vector2 end =
            lineEndpoint(eSegmentEndpointEdge[j], eSegmentEndpointIndex[j]);
        lines.emplace_back(std::make_pair(start, end));
    }

    // add edges from c to d and the opposite vertex
    for (size_t i = 0; i < c; ++i) {
        Vector2 start = z + (double)(i + 1) / ((double)c + 1) * (w - z);
        if (!(cSegmentEndpointEdge[i] == 'd' ||
              cSegmentEndpointEdge[i] == 'x')) {
            continue;
        }
        Vector2 end =
            lineEndpoint(cSegmentEndpointEdge[i], cSegmentEndpointIndex[i]);
        lines.emplace_back(std::make_pair(start, end));
    }

    // add edges from d to the opposite vertex
    for (size_t i = 0; i < d; ++i) {
        Vector2 start = w + (double)(i + 1) / ((double)d + 1) * (x - w);
        if (dSegmentEndpointEdge[i] != 'z') {
            continue;
        }
        Vector2 end = z;
        lines.emplace_back(std::make_pair(start, end));
    }

    size_t numericalIntersections = 0;
    for (std::pair<Vector2, Vector2> line : lines) {
        Vector2 start = line.first;
        Vector2 end   = line.second;
        double iTime =
            intersectionTime(start, end, Vector2{0, -1}, Vector2{0, 1});
        if (iTime > 0 && iTime < 1) {
            numericalIntersections += 1;
        }
    }

    return numericalIntersections;
}

size_t countVerticalEdges(size_t a, size_t b, size_t c, size_t d, size_t e) {
    int topStart = b;
    int topEnd   = positivePart(e - a);

    int bottomStart = c;
    int bottomEnd   = positivePart(e - d);

    int start = std::max(topStart, bottomStart);
    int end   = std::min(topEnd, bottomEnd);

    return positivePart(end - start);
}

void testCoords(size_t a, size_t b, size_t c, size_t d, size_t e) {
    std::array<size_t, 5> neighborhoodNormalCoordinates{e, a, b, c, d};
    size_t f          = flipNormalCoordinate(neighborhoodNormalCoordinates);
    size_t numericalF = computeNumericalEdgeIntersectionNumber(a, b, c, d, e);
    if (e == 0) numericalF += 1;


    if (numericalF != f) {
        // run function in verbose mode
        flipNormalCoordinate(neighborhoodNormalCoordinates);
        cout << "a: " << a << endl
             << "b: " << b << endl
             << "c: " << c << endl
             << "d: " << d << endl
             << "e: " << e << endl
             << "f: " << f << endl
             << "numericalIntersections: " << numericalF << endl;
    }
    ASSERT_EQ(numericalF, f);
}
