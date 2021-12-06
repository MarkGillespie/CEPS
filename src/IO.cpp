#include "IO.h"

namespace CEPS {

std::vector<std::pair<size_t, double>>
readPrescribedConeAngles(std::string filename) {
    return readIndexValuePairs(filename);
}

std::vector<std::pair<size_t, double>>
readPrescribedScaleFactors(std::string filename) {
    return readIndexValuePairs(filename, 0);
}

std::vector<std::pair<size_t, double>>
readIndexValuePairs(std::string filename, double defaultValue) {
    std::ifstream inStream(filename, std::ios::binary);

    return readIndexValuePairs(inStream, defaultValue);
}

std::vector<std::pair<size_t, double>>
readIndexValuePairs(std::istream& in, double defaultValue) {

    std::string line;
    bool done = false;
    std::vector<std::pair<size_t, double>> parsedList;
    while (!done && std::getline(in, line)) {
        std::stringstream ss(line);
        std::string token;

        ss >> token;

        if (token == "#") {
            continue;
        } else {
            size_t index;
            double value;

            // parse token to size_t
            sscanf(token.c_str(), "%zu", &index);

            // try to read a value to go with the index
            if (ss >> value) {
                // read succeeded, nothing to do here
            } else {
                value = defaultValue;
            }
            verbose_assert(!std::isnan(value),
                           "Error, failed to parse line '" + line + "'");

            parsedList.push_back(std::make_pair(index, value));
        }
    }

    return parsedList;
}

void writeMeshWithProjectiveTextureCoords(const SimplePolygonMesh& mesh,
                                          const std::vector<Vector3>& pUV,
                                          std::string filename) {

    std::cout << "Writing mesh to: " << filename << std::endl;

    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("failed to open output file " + filename);
    }
    writeMeshWithProjectiveTextureCoords(mesh, pUV, outFile);
    outFile.close();
}

void writeMeshWithProjectiveTextureCoords(const SimplePolygonMesh& mesh,
                                          const std::vector<Vector3>& pUV,
                                          std::ostream& out) {

    // Write header
    out << "#  vertices: " << mesh.vertexCoordinates.size() << std::endl;
    out << "#     faces: " << mesh.polygons.size() << std::endl;
    out << std::endl;

    // Write vertices
    for (Vector3 p : mesh.vertexCoordinates) {
        out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }


    // Write texture coords
    for (size_t iC = 0; iC < pUV.size(); iC++) {
        out << "vt " << pUV[iC].x << " " << pUV[iC].y << " " << pUV[iC].z
            << std::endl;
    }

    // Write faces
    size_t iC = 0;
    for (const std::vector<size_t>& face : mesh.polygons) {
        out << "f";
        for (size_t ind : face) {
            out << " " << (ind + 1) << "/" << (iC + 1);
            iC++;
        }
        out << std::endl;
    }
}

void writeMeshWithOrdinaryTextureCoords(const SimplePolygonMesh& mesh,
                                        const std::vector<Vector3>& pUV,
                                        std::string filename) {

    std::cout << "Writing mesh to: " << filename << std::endl;

    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("failed to open output file " + filename);
    }
    writeMeshWithOrdinaryTextureCoords(mesh, pUV, outFile);
    outFile.close();
}

void writeMeshWithOrdinaryTextureCoords(const SimplePolygonMesh& mesh,
                                        const std::vector<Vector3>& pUV,
                                        std::ostream& out) {

    // Write header
    out << "#  vertices: " << mesh.vertexCoordinates.size() << std::endl;
    out << "#     faces: " << mesh.polygons.size() << std::endl;
    out << std::endl;

    // Write vertices
    for (Vector3 p : mesh.vertexCoordinates) {
        out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }


    // Write texture coords
    for (size_t iC = 0; iC < pUV.size(); iC++) {
        Vector2 uv{pUV[iC].x / pUV[iC].z, pUV[iC].y / pUV[iC].z};
        out << "vt " << uv.x << " " << uv.y << std::endl;
    }

    // Write faces
    size_t iC = 0;
    for (const std::vector<size_t>& face : mesh.polygons) {
        out << "f";
        for (size_t ind : face) {
            out << " " << (ind + 1) << "/" << (iC + 1);
            iC++;
        }
        out << std::endl;
    }
}

void writeMeshWithProjectiveTextureCoords(
    const SimplePolygonMesh& mesh,
    const std::vector<std::array<double, 4>>& pUV, std::string filename) {

    std::cout << "Writing mesh to: " << filename << std::endl;

    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("failed to open output file " + filename);
    }

    writeMeshWithProjectiveTextureCoords(mesh, pUV, outFile);
    outFile.close();
}

void writeMeshWithProjectiveTextureCoords(
    const SimplePolygonMesh& mesh,
    const std::vector<std::array<double, 4>>& pUV, std::ostream& out) {

    // Write headers
    out << "#  vertices: " << mesh.vertexCoordinates.size() << std::endl;
    out << "#     faces: " << mesh.polygons.size() << std::endl;
    out << std::endl;

    // Write vertices
    for (Vector3 p : mesh.vertexCoordinates) {
        out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }


    // Write texture coords
    for (size_t iC = 0; iC < pUV.size(); iC++) {
        out << "vt " << pUV[iC][0] << " " << pUV[iC][1] << " " << pUV[iC][2]
            << " " << pUV[iC][3] << std::endl;
    }

    // Write faces
    for (const std::vector<size_t>& face : mesh.polygons) {
        out << "f";
        for (size_t ind : face) {
            out << " " << (ind + 1) << "/" << (ind + 1);
        }
        out << std::endl;
    }
}

FaceData<Vector2> readFFieldIntrinsic(std::string filename,
                                      ManifoldSurfaceMesh& mesh,
                                      VertexPositionGeometry& geo) {

    std::vector<double> angles;
    std::ifstream inStream(filename, std::ios::binary);

    std::string line;
    bool done = false;
    while (!done && std::getline(inStream, line)) {
        std::stringstream ss(line);
        std::string token;

        ss >> token;

        if (token == "crossfield_angles") {
            std::getline(inStream, line);
            ss = std::stringstream(line);

            double angle;
            size_t nF = mesh.nFaces();
            for (size_t iF = 0; iF < nF; ++iF) {
                ss >> angle;
                angles.push_back(angle);
            }
            done = true;
        }
    }

    if (!done) {
        std::cerr << "Didn't find crossfield angles?" << std::endl;
    }

    FaceData<Vector2> ffieldData(mesh);
    geo.requireCornerAngles();
    FaceData<size_t> fIdx = mesh.getFaceIndices();
    for (Face f : mesh.faces()) {
        double theta  = geo.cornerAngles[f.halfedge().corner()];
        ffieldData[f] = Vector2::fromAngle(theta + angles[fIdx[f]]);
    }

    return ffieldData;
}

std::vector<std::pair<Vertex, double>>
readFFieldCones(std::string filename, ManifoldSurfaceMesh& mesh,
                VertexPositionGeometry& geo) {
    FaceData<Vector2> crossField = readFFieldIntrinsic(filename, mesh, geo);
    VertexData<int> vertexIndices =
        computeVertexIndex(mesh, geo, crossField, 4);

    Vertex vBad;
    std::vector<std::pair<Vertex, double>> cones;
    for (Vertex v : mesh.vertices()) {
        if (vertexIndices[v] != 0 && !v.isBoundary()) {
            cones.push_back(std::make_pair(v, M_PI / 2. * vertexIndices[v]));
        }
        if (vertexIndices[v] > 4) vBad = v;
    }

    return cones;
}
} // namespace CEPS
