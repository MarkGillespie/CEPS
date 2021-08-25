#include "Utils.h"

verbose_runtime_error::verbose_runtime_error(const std::string& arg,
                                             const char* file, int line)
    : std::runtime_error(arg) {

    std::string filePath(file);

    // stolen from polyscope/utilities.cpp
    size_t startInd = 0;
    for (std::string sep : {"/", "\\"}) {
        size_t pos = filePath.rfind(sep);
        if (pos != std::string::npos) {
            startInd = std::max(startInd, pos + 1);
        }
    }

    std::string niceName = filePath.substr(startInd, filePath.size());

    std::ostringstream o;
    o << arg << " At " << niceName << ":" << line;
    msg = o.str();
}

namespace CEPS {
Vertex src(Edge e) { return e.halfedge().vertex(); }
Vertex dst(Edge e) { return e.halfedge().next().vertex(); }

int positivePart(int x) { return fmax(x, 0); }
int negativePart(int x) { return fmin(x, 0); }

Vector2 bary2(double t) {
    if (t <= 0) {
        return Vector2{0, 1};
    } else if (t >= 1) {
        return Vector2{1, 0};
    } else {
        return Vector2{t, 1 - t};
    }
}

Vector2 normalizeBary(Vector2 b) { return b / (b.x + b.y); }
Vector3 normalizeBary(Vector3 b) { return b / (b.x + b.y + b.z); }

// Returns the intersection of the line a1-a2 with line b1-b2 in barycentric
// coordinates along line a1-a2 (i.2. t*a1 + (1-t)a2 lies on line b1-b2 as
// well)
double intersectionTime(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2) {
    // TODO: are these checks time-consuming?
    if ((a1 - b1).norm() < 1e-12 || (a1 - b2).norm() < 1e-12) {
        return 1;
    } else if ((a2 - b1).norm() < 1e-12 || (a2 - b2).norm() < 1e-12) {
        return 0;
    } else if ((a1 - a2).norm() < 1e-12 && (b1 - b2).norm() < 1e-12) {
        // If a and b coincided, we would have caught that in a previous case.
        // So the "lines" don't intersect, and we just return a big number
        return 50;
    } else if ((b1 - b2).norm() < 1e-12) {
        Vector2 dA = (a2 - a1).normalize();
        return 1 - dot(b1 - a1, dA) / dot(a2 - a1, dA);
    } else if ((a1 - a2).norm() < 1e-12) {
        Vector2 dB = (b2 - b1).normalize();
        double t   = dot(a1 - b1, dB) / dot(b2 - b1, dB);

        if (0 <= t && t <= 1) {
            Vector2 projAB = b1 + t * dB;
            if ((projAB - a1).norm() < 1e-12) {
                return 0.5;
            } else {
                return 50;
            }
        } else {
            return 50;
        }

    } else {
        double m1 = a2.x - a1.x;
        double m2 = b2.x - b1.x;
        double m3 = a2.y - a1.y;
        double m4 = b2.y - b1.y;
        double vx = b1.x - a1.x;
        double vy = b1.y - a1.y;
        double t  = (m2 * vy - m4 * vx) / (m2 * m3 - m1 * m4);

        return 1 - t;
    }
}

// Returns the barycentric coordinate for pt along line a-b
// i.e. minimizes \|b + (a-b)t - pt\|_2^2
// i.e pt = ta + (1-t)b
// Clamps t to the interval [0, 1]
double barycentricCoordinate(Vector2 pt, Vector2 a, Vector2 b) {
    double t = dot(b - a, b - pt) / norm2(b - a);
    if (t >= 0 && t <= 1) {
        return t;
    } else {
        if (t < -1e-12 || t > 1 + 1e-12) {
            // cerr << "error Barycentric coordinate " << t
            //      << " invalid [barcentricCoordinate(Vector2 ...)]" <<
            //      endl;
        }
        if (t < 0) {
            return 0;
        } else {
            return 1;
        }
    }
}

double lobachevsky(double a) { return 0.5 * clausen(2 * a); }
double safeSqrt(double x) { return sqrt(fmax(x, 0)); }

void symmetrizeVizRange(polyscope::SurfaceVertexIsolatedScalarQuantity& q) {
    q.vizRangeLow  = fmin(q.dataRangeLow, -q.dataRangeHigh);
    q.vizRangeHigh = fmax(q.dataRangeHigh, -q.dataRangeLow);
}
} // namespace CEPS


string getFilename(string filePath, bool withExtension, char seperator) {
    // Get last dot position
    size_t dotPos = filePath.rfind('.');
    size_t sepPos = filePath.rfind(seperator);

    if (sepPos != string::npos) {
        return filePath.substr(
            sepPos + 1,
            filePath.size() -
                (withExtension || dotPos != std::string::npos ? 1 : dotPos));
    }
    return "";
}
