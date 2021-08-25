#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/SurfaceProjectiveParameterizationQuantity.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "ConeFlattening.h"
#include "IO.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace CEPS;

//== Input data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::vector<size_t> imaginaryFaceIndices;
std::vector<std::pair<Vertex, double>> prescribedCurvatures;
std::vector<std::pair<Vertex, double>> prescribedScaleFactors;

//== Result data
ParameterizationResult result;

//== Visualization data
polyscope::SurfaceMesh* psMesh;

const int IM_STR_LEN = 128;
static char meshSaveName[IM_STR_LEN];

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Place greedy cones")) {
        prescribedScaleFactors.clear();
        prescribedCurvatures.clear();

        std::vector<Vertex> coneVertices = placeCones(*mesh, *geometry);

        // list of cone indices, with a dummy second variable since polyscope
        // wants to render a function on the cones
        std::vector<std::pair<size_t, int>> coneIndices;
        for (Vertex v : coneVertices) {
            coneIndices.emplace_back(v.getIndex(), 0);
            prescribedScaleFactors.emplace_back(v, 0);
        }
        psMesh->addVertexCountQuantity("cones", coneIndices);
    }

    if (ImGui::Button("Uniformize with greedy cones")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterizeWithGreedyCones(*mesh, *geometry, viz,
                                             checkInjectivity);
        psMesh->setEnabled(false);
    }

    if (ImGui::Button("Uniformize")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterize(*mesh, *geometry, prescribedScaleFactors,
                              prescribedCurvatures, imaginaryFaceIndices, viz,
                              checkInjectivity);
        psMesh->setEnabled(false);
    }

    ImGui::Separator();
    ImGui::InputText("###PtexturedMeshSaveName", meshSaveName,
                     IM_ARRAYSIZE(meshSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save textured mesh")) {
        writeMeshWithProjectiveTextureCoords(
            result.mesh, result.param, std::string(meshSaveName) + ".obj");
        saveMatrix(result.interpolationMatrix,
                   std::string(meshSaveName) + ".spmat");
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Geometry program");
    args::Positional<std::string> meshFilename(parser, "mesh",
                                               "Mesh to be processed.");
    args::ValueFlag<std::string> curvaturesFilename(
        parser, "string", "curvatures filename", {"curvatures"});
    args::ValueFlag<std::string> scaleFactorsFilename(
        parser, "string", "scaleFactors filename", {"scaleFactors"});
    args::ValueFlag<std::string> ffieldFilename(parser, "string",
                                                "ffield filename", {"ffield"});
    args::Flag useExactCones(parser, "exactCones",
                             "Use exact cones from ffield (no lumping)",
                             {"exactCones"});
    args::Flag minimalAreaDistortion(parser, "minimalAreaDistortion",
                                     "Apply minimal-area-distortion boundary "
                                     "conditions when others are not specified",
                                     {"minimalAreaDistortion"});

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help& h) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    if (!meshFilename) {
        std::cout << "Please provide a mesh file as input" << std::endl;
        return 1;
    }

    std::string filename = args::get(meshFilename);

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(filename);

    std::string nicename = polyscope::guessNiceNameFromPath(filename);

    std::string saveNameGuess = nicename + "_ceps";
    // Initialize ImGui's C-string to our niceName
    // https://stackoverflow.com/a/347959
    // truncate saveNameGuess to fit in ImGui's string
    if (saveNameGuess.length() + 1 > IM_STR_LEN)
        saveNameGuess.resize(IM_STR_LEN - 1);
    // copy over string contents
    std::copy(saveNameGuess.begin(), saveNameGuess.end(), meshSaveName);
    // null-terminate string
    meshSaveName[saveNameGuess.size()] = '\0';

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        "input_mesh", geometry->inputVertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));

    // Read prescribed curvatures and scale factors if present
    if (curvaturesFilename) {
        std::vector<std::pair<size_t, double>> prescribedCurvatureInput =
            readPrescribedConeAngles(args::get(curvaturesFilename));
        auto q = psMesh->addVertexIsolatedScalarQuantity(
            "prescribed curvatures", prescribedCurvatureInput);
        symmetrizeVizRange(*q);
        q->setEnabled(true);
        for (const std::pair<size_t, double>& c : prescribedCurvatureInput)
            prescribedCurvatures.emplace_back(mesh->vertex(c.first), c.second);
    }
    if (scaleFactorsFilename) {
        std::vector<std::pair<size_t, double>> prescribedScaleFactorInput =
            readPrescribedScaleFactors(args::get(scaleFactorsFilename));
        auto q = psMesh->addVertexIsolatedScalarQuantity(
            "prescribed scaleFactors", prescribedScaleFactorInput);
        symmetrizeVizRange(*q);
        q->setEnabled(true);
        for (const std::pair<size_t, double>& c : prescribedScaleFactorInput)
            prescribedScaleFactors.emplace_back(mesh->vertex(c.first),
                                                c.second);
    }
    if (ffieldFilename) {

        std::vector<std::pair<Vertex, double>> ffieldCones =
            readFFieldCones(args::get(ffieldFilename), *mesh, *geometry);

        auto processedCones = args::get(useExactCones)
                                  ? ffieldCones
                                  : lumpCones(*mesh, ffieldCones);

        std::vector<std::pair<size_t, double>> coneIndices;
        for (const std::pair<Vertex, double>& c : processedCones)
            coneIndices.emplace_back(c.first.getIndex(), c.second);
        auto q = psMesh->addVertexIsolatedScalarQuantity("ffield cones",
                                                         coneIndices);
        symmetrizeVizRange(*q);
        q->setEnabled(true);

        prescribedCurvatures.insert(std::end(prescribedCurvatures),
                                    std::begin(processedCones),
                                    std::end(processedCones));
    }

    if (args::get(minimalAreaDistortion)) {
        // Identify which vertices are already fixed
        VertexData<bool> prescribed(*mesh);
        for (const std::pair<Vertex, double> cone : prescribedCurvatures)
            prescribed[cone.first] = true;
        for (const std::pair<Vertex, double> cone : prescribedScaleFactors)
            prescribed[cone.first] = true;

        // Set u=0 at all other boundary vertices
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                if (!prescribed[v]) prescribedScaleFactors.emplace_back(v, 0);
            }
        }
    }

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
