#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "IO.h"
#include "Logger.h"
#include "SphericalUniformization.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace CEPS;

//== Input data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

//== Result data
SphericalUniformizationResult result;

//== Visualization data
polyscope::SurfaceMesh* psMesh;

const int IM_STR_LEN = 128;
static char meshSaveName[IM_STR_LEN];

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Uniformize")) {
        bool viz = true;
        result   = sphericalUniformize(*mesh, *geometry, Vertex(), viz);
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
    args::ArgumentParser parser("spherical uniformization demo");
    args::Positional<std::string> meshFilename(parser, "mesh",
                                               "Mesh to be processed.");
    args::ValueFlag<std::string> outputMeshFilename(
        parser, "string", "output mesh filename", {"outputMeshFilename"});
    args::ValueFlag<std::string> outputMatrixFilename(
        parser, "string", "output matrix filename", {"outputMatrixFilename"});
    args::ValueFlag<std::string> outputLogFilename(
        parser, "string", "output log filename", {"outputLogFilename"});
    args::Flag viz(parser, "viz", "Use polyscope GUI", {"viz"});

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

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(filename);

    // center mesh
    Vector3 mean = Vector3::zero();
    for (Vertex v : mesh->vertices()) mean += geometry->inputVertexPositions[v];
    mean /= mesh->nVertices();
    for (Vertex v : mesh->vertices()) geometry->inputVertexPositions[v] -= mean;

    std::string nicename = polyscope::guessNiceNameFromPath(filename);

    std::string saveNameGuess = nicename + "_spherical";
    // Initialize ImGui's C-string to our niceName
    // https://stackoverflow.com/a/347959
    // truncate saveNameGuess to fit in ImGui's string
    if (saveNameGuess.length() + 1 > IM_STR_LEN)
        saveNameGuess.resize(IM_STR_LEN - 1);
    // copy over string contents
    std::copy(saveNameGuess.begin(), saveNameGuess.end(), meshSaveName);
    // null-terminate string
    meshSaveName[saveNameGuess.size()] = '\0';

    if (args::get(viz)) {
        // Initialize polyscope
        polyscope::init();

        // Set the callback function
        polyscope::state::userCallback = myCallback;

        // Register the mesh with polyscope
        psMesh = polyscope::registerSurfaceMesh(
            "input_mesh", geometry->inputVertexPositions,
            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        std::cout << "Loaded mesh " << filename << std::endl;
        std::cout << "nBoundaryLoops: " << mesh->nBoundaryLoops() << std::endl;
        std::cout << "Genus: " << mesh->genus() << std::endl;
        std::cout << "Euler characteristic: " << mesh->eulerCharacteristic()
                  << std::endl;
        std::cout << "Connected components: " << mesh->nConnectedComponents()
                  << std::endl;

        // Give control to the polyscope gui
        polyscope::show();
    } else {
        if (outputLogFilename) {
            std::ofstream out;

            // std::ios::trunc ensures that we overwrite old versions
            out.open(args::get(outputLogFilename), std::ios::trunc);
            if (out.is_open()) {
                // output a temporary symbol so we can tell if the program
                // crashes before writing the real log
                out << ":'(" << std::endl;
                out.close();
            } else {
                std::cout << "Error: failed to write to "
                          << args::get(outputLogFilename) << vendl;
            }
        }

        Logger logger;
        bool logStats = false;
        if (outputLogFilename) {
            logStats = true;
            logger.log("name", nicename);
            logger.log("nVertices", mesh->nVertices());
        }

        double duration;
        std::clock_t start = std::clock();

        bool viz = false;
        result   = sphericalUniformize(*mesh, *geometry, Vertex(), viz);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;


        if (outputMeshFilename) {
            writeMeshWithProjectiveTextureCoords(result.mesh, result.param,
                                                 args::get(outputMeshFilename));
        }

        if (outputMatrixFilename) {
            saveMatrix(result.interpolationMatrix,
                       args::get(outputMatrixFilename));
        }

        if (logStats) {
            logger.log("duration", duration);
            logger.writeLog(args::get(outputLogFilename));
        }
    }

    return EXIT_SUCCESS;
}
