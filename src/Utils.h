#pragma once

#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_count_quantity.h"

#include "Clausen.h"

#include <iostream>
#include <sstream>
#include <string>

using std::string;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Custom runtime error which will report the offending file and line number
class verbose_runtime_error : public std::runtime_error {
    std::string msg;

  public:
    verbose_runtime_error(const std::string& arg, const char* file, int line);
    ~verbose_runtime_error() throw() {}
    const char* what() const throw() { return msg.c_str(); }
};

#define throw_verbose_runtime_error(arg)                                       \
    throw verbose_runtime_error(arg, __FILE__, __LINE__);

#define verbose_assert(arg, msg)                                               \
    if (!(arg)) throw_verbose_runtime_error(msg);


namespace CEPS {
enum class GeometryType { EUCLIDEAN, HYPERBOLIC };
Vertex src(Edge e);
Vertex dst(Edge e);

int positivePart(int x);
int negativePart(int x);

Vector2 bary2(double t);

Vector2 normalizeBary(Vector2 b);
Vector3 normalizeBary(Vector3 b);

// Returns the intersection of the line a1-a2 with line b1-b2 in barycentric
// coordinates along line a1-a2 (i.2. t*a1 + (1-t)a2 lies on line b1-b2 as
// well)
double intersectionTime(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2);

// Returns the barycentric coordinate for pt along line a-b
// i.e. minimizes \|b + (a-b)t - pt\|_2^2
// i.e pt = ta + (1-t)b
// Clamps t to the interval [0, 1]
double barycentricCoordinate(Vector2 pt, Vector2 a, Vector2 b);

double lobachevsky(double a);
double safeSqrt(double x);

void symmetrizeVizRange(polyscope::SurfaceVertexIsolatedScalarQuantity& q);

} // namespace CEPS

//== Convenient printing helpers

string getFilename(string filePath, bool withExtension = true,
                   char seperator = '/');

// Verbose endl
#define vendl                                                                  \
    "\t\t(" << getFilename(__FILE__) << ":" << __LINE__ << ")" << std::endl

// == Horrible macro to print variable name + variable
// https://stackoverflow.com/a/6623090
#define WATCHSTR_WNAME(os, name, a)                                            \
    do {                                                                       \
        (os) << (name) << " is value " << (a) << vendl;                        \
    } while (false)

#define WATCHSTR(os, a) WATCHSTR_WNAME((os), #a, (a))
#define WATCH(a) WATCHSTR_WNAME(std::cout, #a, (a))
#define WATCH2(a, b)                                                           \
    WATCH(a);                                                                  \
    WATCH(b)
#define WATCH3(a, b, c)                                                        \
    WATCH(a);                                                                  \
    WATCH2(b, c)
#define WATCH4(a, b, c, d)                                                     \
    WATCH(a);                                                                  \
    WATCH3(b, c, d)
#define WATCH5(a, b, c, d, e)                                                  \
    WATCH(a);                                                                  \
    WATCH4(b, c, d, e)

#define HERE()                                                                 \
    do {                                                                       \
        std::cout << "HERE at " << __LINE__ << " in " << getFilename(__FILE__) \
                  << ":" << __FUNCTION__ << std::endl;                         \
    } while (false)
#define HERE1(name)                                                            \
    do {                                                                       \
        std::cout << "HERE (" << (name) << ") at " << __LINE__ << " in "       \
                  << getFilename(__FILE__) << ":" << __FUNCTION__              \
                  << std::endl;                                                \
    } while (false)

namespace CEPS {
namespace ImplementationDetails {

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data);

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& data);

template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& data);

} // namespace ImplementationDetails
} // namespace CEPS


#include "Utils.ipp"
