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
bool verbose = false;

//== Result data
SphericalUniformizationResult result;

//== Visualization data
polyscope::SurfaceMesh* psMesh;

const int IM_STR_LEN = 128;
static char meshSaveName[IM_STR_LEN];

enum class InterpolationType { PROJECTIVE, LINEAR };
int currentInterpolationType                       = 0;
const char* prettyInterpolationTypeOptions[]       = {"Homogeneous", "Linear"};
const InterpolationType interpolationTypeOptions[] = {
    InterpolationType::PROJECTIVE, InterpolationType::LINEAR};

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Uniformize")) {
        bool viz = true;
        result = sphericalUniformize(*mesh, *geometry, Vertex(), viz, verbose);
        psMesh->setEnabled(false);
    }

    ImGui::Separator();
    ImGui::InputText("###PtexturedMeshSaveName", meshSaveName,
                     IM_ARRAYSIZE(meshSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save textured mesh")) {
        switch (interpolationTypeOptions[currentInterpolationType]) {
        case InterpolationType::PROJECTIVE:
            writeMeshWithProjectiveTextureCoords(
                result.mesh, result.param, std::string(meshSaveName) + ".obj");
            break;
        case InterpolationType::LINEAR:
            writeMeshWithOrdinaryTextureCoords(
                result.mesh, result.param, std::string(meshSaveName) + ".obj");
            break;
        }
        saveMatrix(result.interpolationMatrix,
                   std::string(meshSaveName) + ".spmat");
    }
    ImGui::Combo("Saved Texture Type", &currentInterpolationType,
                 prettyInterpolationTypeOptions,
                 IM_ARRAYSIZE(interpolationTypeOptions));
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("spherical uniformization demo");
    args::Positional<std::string> meshFilename(
        parser, "mesh", "Mesh to be processed (required).");
    args::ValueFlag<std::string> outputMeshFilename(
        parser, "string",
        "file to save output mesh to, along with homogeneous texture "
        "coordinates",
        {"outputMeshFilename"});
    args::ValueFlag<std::string> outputLinearTextureFilename(
        parser, "string",
        "file to save output mesh to, along with linear texture coordinates "
        "(aka ordinary uv coordinates)",
        {"outputLinearTextureFilename"});
    args::ValueFlag<std::string> outputMatrixFilename(
        parser, "string", "file to save the output interpolation matrix to",
        {"outputMatrixFilename"});
    args::ValueFlag<std::string> outputLogFilename(
        parser, "string", "file to save logs to", {"outputLogFilename"});
    args::Flag viz(parser, "viz", "Use polyscope GUI", {"viz"});
    args::ValueFlag<std::string> beVerbose(
        parser, "verbose",
        "[y/n] Print out progress information (default "
        "true in GUI mode, false otherwise)",
        {"verbose"});
    args::Flag version(parser, "version", "Display version number",
                       {'v', "version"});
    args::HelpFlag help(parser, "help", "Display this help menu",
                        {'h', "help"});

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

    if (version) {
        std::cout << "spherical_uniformize version 1.2" << std::endl;
        return 0;
    }

    if (!meshFilename) {
        std::cout << "Please provide a mesh file as input" << std::endl;
        std::cout << parser;
        return 1;
    }

    verbose = args::get(viz);
    if (beVerbose) {
        std::string verbosity = args::get(beVerbose);
        std::transform(verbosity.begin(), verbosity.end(), verbosity.begin(),
                       ::toupper);
        if (verbosity == "Y" || verbosity == "TRUE") {
            verbose = true;
        } else if (verbosity == "N" || verbosity == "FALSE") {
            verbose = false;
        } else {
            std::cout << "Unrecognized verbosity option '" << verbosity << "'"
                      << std::endl;
            std::cout << "\t please use true/false or y/n" << std::endl;
            std::cout << parser;
            return 1;
        }
    }

    std::string filename = args::get(meshFilename);

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(filename);
    verbose_assert(mesh->nConnectedComponents() == 1, "mesh must be connected");
    verbose_assert(mesh->genus() == 0, "mesh must have genus 0 (but it is " +
                                           std::to_string(mesh->genus()) + ")");
    verbose_assert(mesh->nBoundaryLoops() == 0, "mesh must not have boundary");

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
        result = sphericalUniformize(*mesh, *geometry, Vertex(), viz, verbose);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;


        if (outputMeshFilename) {
            writeMeshWithProjectiveTextureCoords(result.mesh, result.param,
                                                 args::get(outputMeshFilename));
        }

        if (outputLinearTextureFilename) {
            writeMeshWithOrdinaryTextureCoords(
                result.mesh, result.param,
                args::get(outputLinearTextureFilename));
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
