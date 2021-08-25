#include "geometrycentral/surface/simple_polygon_mesh.h"

#include "GeometryUtils.h"
#include "Utils.h"

namespace CEPS {

std::vector<std::pair<size_t, double>>
readPrescribedConeAngles(std::string filename);

std::vector<std::pair<size_t, double>>
readPrescribedScaleFactors(std::string filename);

//== Read in a list of index-value pairs from a file.
// Each pair must appear on its own line.
// Lines beginning with '#' are ignored.
// If defaultValue is set to something other than NaN, then indices which do not
// have values are assigned the default value. Otherwise, throws an error if an
// index appears without a corresponding value.
std::vector<std::pair<size_t, double>> readIndexValuePairs(
    std::string filename,
    double defaultValue = std::numeric_limits<double>::quiet_NaN());

std::vector<std::pair<size_t, double>> readIndexValuePairs(
    std::istream& in,
    double defaultValue = std::numeric_limits<double>::quiet_NaN());

//== Save mesh to obj with projective texture coordinates
// pUV gives UV coordinates at *corners*
void writeMeshWithProjectiveTextureCoords(const SimplePolygonMesh& mesh,
                                          const std::vector<Vector3>& pUV,
                                          std::string filename);

void writeMeshWithProjectiveTextureCoords(const SimplePolygonMesh& mesh,
                                          const std::vector<Vector3>& pUV,
                                          std::ostream& out);

// pUV gives UV coordinates at *vertices*.
void writeMeshWithProjectiveTextureCoords(
    const SimplePolygonMesh& mesh,
    const std::vector<std::array<double, 4>>& pUV, std::string filename);

void writeMeshWithProjectiveTextureCoords(
    const SimplePolygonMesh& mesh,
    const std::vector<std::array<double, 4>>& pUV, std::ostream& out);

//== Save a sparse matrix as a list of triplets. Uses 1-indexing
template <typename T>
void saveMatrix(Eigen::SparseMatrix<T>& matrix, std::string filename);

template <typename T>
void saveMatrix(Eigen::SparseMatrix<T>& matrix, std::ostream& out);


//== Read cross fields from MPZ
// Read the intrinsic cross field description. Output cross field is encoded by
// an arbitrary choice of representative vector per face (i.e. not raised to the
// 4th power)
FaceData<Vector2> readFFieldIntrinsic(std::string filename,
                                      ManifoldSurfaceMesh& mesh,
                                      VertexPositionGeometry& geo);

// Extract the cone angles (i.e. singularities) from an MPZ cross field
std::vector<std::pair<Vertex, double>>
readFFieldCones(std::string filename, ManifoldSurfaceMesh& mesh,
                VertexPositionGeometry& geo);
} // namespace CEPS

#include "IO.ipp"
