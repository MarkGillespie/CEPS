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
#include "Logger.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace CEPS;

//== Input data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::vector<size_t> imaginaryFaceIndices;
std::vector<std::pair<Vertex, double>> prescribedCurvatures;
std::vector<std::pair<Vertex, double>> prescribedScaleFactors;
double uTol  = 5; // maximum allowed area distortion in greedy cone placement
bool verbose = false;

//== Result data
ParameterizationResult result;

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
    if (ImGui::Button("Place greedy cones")) {
        prescribedScaleFactors.clear();
        prescribedCurvatures.clear();

        std::vector<Vertex> coneVertices =
            placeCones(*mesh, *geometry, uTol, 4, 4, verbose);

        // list of cone indices, with a dummy second variable since polyscope
        // wants to render a function on the cones
        std::vector<std::pair<size_t, int>> coneIndices;
        for (Vertex v : coneVertices) {
            coneIndices.emplace_back(v.getIndex(), 0);
            prescribedScaleFactors.emplace_back(v, 0);
        }

        // Minimal area distortion boundary conditions.
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                prescribedScaleFactors.emplace_back(v, 0);
            }
        }

        psMesh->addVertexCountQuantity("cones", coneIndices);
    }

    if (ImGui::Button("Uniformize with greedy cones")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterizeWithGreedyCones(*mesh, *geometry, viz,
                                             checkInjectivity, uTol, verbose);
        psMesh->setEnabled(false);

        std::cout << "nInvertedTriangles: " << result.nFlippedTriangles
                  << "\tnZeroAreaTriangles: " << result.nZeroAreaTriangles
                  << std::endl;
    }

    if (ImGui::Button("Uniformize")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterize(*mesh, *geometry, prescribedScaleFactors,
                              prescribedCurvatures, imaginaryFaceIndices, viz,
                              checkInjectivity, verbose);
        psMesh->setEnabled(false);

        std::cout << "nInvertedTriangles: " << result.nFlippedTriangles
                  << "\tnZeroAreaTriangles: " << result.nZeroAreaTriangles
                  << std::endl;
    }

    ImGui::Separator();
    ImGui::InputText("###PtexturedMeshSaveName", meshSaveName,
                     IM_ARRAYSIZE(meshSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save Texture")) {
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
    args::ArgumentParser parser("Conformal cone flattener");
    args::Positional<std::string> meshFilename(
        parser, "mesh", "Mesh to be processed (required).");
    args::ValueFlag<std::string> curvaturesFilename(
        parser, "string", "curvatures filename", {"curvatures"});
    args::ValueFlag<std::string> scaleFactorsFilename(
        parser, "string", "scaleFactors filename", {"scaleFactors"});
    args::ValueFlag<std::string> ffieldFilename(parser, "string",
                                                "ffield filename", {"ffield"});
    args::ValueFlag<double> greedyConeMaxU(
        parser, "double", "max allowed scale factor in greedy cone placement",
        {"greedyConeMaxU"});

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
        parser, "string", "file to save logs (timing + injectivity) to",
        {"outputLogFilename"});

    args::Flag useExactCones(parser, "exactCones",
                             "Use exact cones from ffield (no lumping)",
                             {"exactCones"});
    args::Flag noFreeBoundary(
        parser, "noFreeBoundary",
        "Do not apply minimal-area-distortion boundary "
        "conditions when others are not specified. Useful when prescribing, "
        "e.g., polygonal boundary conditions",
        {"noFreeBoundary"});
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
        std::cout << "parameterize version 1.1" << std::endl;
        return 0;
    }

    if (!meshFilename) {
        std::cout << "Please provide a mesh file as input." << std::endl;
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

    bool loadedCones = false;
    // Read prescribed curvatures and scale factors if present
    std::vector<std::pair<size_t, double>> prescribedCurvatureInput,
        prescribedScaleFactorInput, ffieldConeInput;

    if (curvaturesFilename) {
        loadedCones = true;
        prescribedCurvatureInput =
            readPrescribedConeAngles(args::get(curvaturesFilename));
        for (const std::pair<size_t, double>& c : prescribedCurvatureInput)
            prescribedCurvatures.emplace_back(mesh->vertex(c.first), c.second);
    }
    if (scaleFactorsFilename) {
        loadedCones = true;
        prescribedScaleFactorInput =
            readPrescribedScaleFactors(args::get(scaleFactorsFilename));
        for (const std::pair<size_t, double>& c : prescribedScaleFactorInput)
            prescribedScaleFactors.emplace_back(mesh->vertex(c.first),
                                                c.second);
    }
    if (ffieldFilename) {
        loadedCones = true;
        std::vector<std::pair<Vertex, double>> ffieldCones =
            readFFieldCones(args::get(ffieldFilename), *mesh, *geometry);

        auto processedCones = args::get(useExactCones)
                                  ? ffieldCones
                                  : lumpCones(*mesh, ffieldCones);

        for (const std::pair<Vertex, double>& c : processedCones)
            ffieldConeInput.emplace_back(c.first.getIndex(), c.second);

        prescribedCurvatures.insert(std::end(prescribedCurvatures),
                                    std::begin(processedCones),
                                    std::end(processedCones));
    }

    // Set u=0 at any remaining boundary vertices unless instructed not to
    if (!args::get(noFreeBoundary)) {
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

    if (greedyConeMaxU) {
        uTol = args::get(greedyConeMaxU);
    }

    if (args::get(viz)) {
        // Initialize polyscope
        polyscope::init();

        // Set the callback function
        polyscope::state::userCallback = myCallback;

        // Register the mesh with polyscope
        psMesh = polyscope::registerSurfaceMesh(
            "input_mesh", geometry->inputVertexPositions,
            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        if (viz) {
            if (!prescribedCurvatureInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "prescribed curvatures", prescribedCurvatureInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }

            if (!prescribedScaleFactorInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "prescribed scaleFactors", prescribedScaleFactorInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }

            if (!ffieldConeInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "ffield cones", ffieldConeInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }
        }

        std::cout << "Loaded mesh " << filename << std::endl;
        std::cout << "nBoundaryLoops: " << mesh->nBoundaryLoops() << std::endl;
        std::cout << "Genus: " << mesh->genus() << std::endl;
        std::cout << "Euler characteristic: " << mesh->eulerCharacteristic()
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
        if (outputLogFilename) logStats = true;

        if (logStats) {
            logger.log("name", nicename);
            logger.log("nVertices", mesh->nVertices());
        }

        double duration;
        if (loadedCones) {
            std::clock_t start = std::clock();

            bool viz              = false;
            bool checkInjectivity = logStats;
            result = parameterize(*mesh, *geometry, prescribedScaleFactors,
                                  prescribedCurvatures, imaginaryFaceIndices,
                                  viz, checkInjectivity, verbose);

            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        } else {
            std::clock_t start = std::clock();

            bool viz                  = false;
            bool checkInjectivity     = logStats;
            std::vector<Vertex> cones = placeCones(*mesh, *geometry, uTol);
            result = parameterizeWithGivenCones(*mesh, *geometry, cones, viz,
                                                checkInjectivity, verbose);

            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        }

        if (logStats) {
            logger.log("duration", duration);
            logger.log("nFlippedTriangles", result.nFlippedTriangles);
            logger.log("nZeroAreaTriangles", result.nZeroAreaTriangles);
        }

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

        if (outputLogFilename) {
            logger.writeLog(args::get(outputLogFilename));
        }
    }

    return EXIT_SUCCESS;
}
