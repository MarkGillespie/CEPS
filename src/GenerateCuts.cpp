#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

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
std::vector<Vertex> cones;
std::set<Edge> cuts;
double uTol  = 5;  // maximum allowed area distortion in greedy cone placement
int nCones   = 10; // number of random cones to place
bool verbose = false;

//== Output data
std::unique_ptr<ManifoldSurfaceMesh> cutMesh;
std::unique_ptr<VertexPositionGeometry> cutGeometry;

//== Visualization data
polyscope::SurfaceMesh* psMesh;

const int IM_STR_LEN = 128;
static char meshSaveName[IM_STR_LEN];

void psDrawCones() {
    // list of cone indices, with a dummy second variable since polyscope
    // wants to render a function on the cones
    std::vector<std::pair<size_t, int>> coneIndices;
    for (Vertex v : cones) {
        coneIndices.emplace_back(v.getIndex(), 0);
    }

    psMesh->addVertexCountQuantity("cones", coneIndices);
}

void psDrawCuts() {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeInds;

    nodePositions.reserve(cuts.size() * 2);
    edgeInds.reserve(cuts.size());
    geometry->requireVertexPositions();
    for (Edge ij : cuts) {
        // if (ij == Edge()) {
        //     std::cout << "WARNING: invalid edge" << std::endl;
        //     continue;
        // }
        Vertex i = ij.halfedge().tailVertex();
        Vertex j = ij.halfedge().tipVertex();
        // WATCH3(i, j, mesh->nVertices());
        // WATCH(geometry->vertexPositions[i]);
        // WATCH(geometry->vertexPositions[j]);
        edgeInds.push_back({nodePositions.size(), nodePositions.size() + 1});
        nodePositions.push_back(geometry->vertexPositions[i]);
        nodePositions.push_back(geometry->vertexPositions[j]);
    }
    geometry->unrequireVertexPositions();

    psMesh->addSurfaceGraphQuantity("cuts", nodePositions, edgeInds);
}

void psDrawCutMesh() {
    polyscope::registerSurfaceMesh(
        "Cut Mesh", cutGeometry->inputVertexPositions,
        cutMesh->getFaceVertexList(), polyscopePermutations(*cutMesh));
}

void generateRandomCones() {
    std::set<size_t> coneIndexSet;
    while (coneIndexSet.size() < std::min(nCones, (int)mesh->nVertices())) {
        coneIndexSet.insert(randomInt(0, mesh->nVertices() - 1));
    }

    cones.clear();
    cones.reserve(coneIndexSet.size());
    for (size_t iV : coneIndexSet) {
        cones.push_back(mesh->vertex(iV));
    }
}

void generateGreedyCones() {
    cones = placeCones(*mesh, *geometry, uTol, 4, 4, verbose);
}

void generateGoodCuts() {
    EdgeData<double> edgeCost(*mesh, 1);
    cuts = ImplementationDetails::goodCuts(*mesh, cones, edgeCost);
}

void generateHomologyGenerators() {
    EdgeData<double> edgeCost(*mesh, 1);
    std::vector<std::vector<Edge>> homologyGenerators =
        ImplementationDetails::computeHomologyGenerators(*mesh, edgeCost);

    cuts.clear();
    for (const std::vector<Edge>& generator : homologyGenerators) {
        for (Edge e : generator) cuts.insert(e);
    }
}

void cutAlongCuts() {
    CornerData<size_t> cutVertexIndices;
    size_t nV;
    std::tie(cutVertexIndices, nV) =
        ImplementationDetails::indexCorners(*mesh, cuts);

    std::vector<Vector3> cutVertexPositions;
    cutVertexPositions.resize(nV);
    std::vector<std::vector<size_t>> cutFaces;
    cutFaces.reserve(mesh->nFaces());

    geometry->requireVertexPositions();
    for (Face f : mesh->faces()) {
        cutFaces.push_back({});
        cutFaces.back().reserve(f.degree());
        for (Corner c : f.adjacentCorners()) {
            cutVertexPositions[cutVertexIndices[c]] =
                geometry->vertexPositions[c.vertex()];
            cutFaces.back().push_back(cutVertexIndices[c]);
        }
    }
    geometry->unrequireVertexPositions();

    std::tie(cutMesh, cutGeometry) =
        makeManifoldSurfaceMeshAndGeometry(cutFaces, cutVertexPositions);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Place greedy cones")) {
        generateGreedyCones();
        psDrawCones();
    }

    ImGui::SliderInt("# Random Cones", &nCones, 1, 100);
    if (ImGui::Button("Place random cones")) {
        generateRandomCones();
        psDrawCones();
    }

    if (ImGui::Button("Generate Cuts")) {
        generateGoodCuts();
        psDrawCuts();
    }
    if (ImGui::Button("Generate Homology Generators")) {
        generateHomologyGenerators();
        psDrawCuts();
    }

    if (ImGui::Button("Cut Mesh Along Cuts")) {
        cutAlongCuts();
        psDrawCutMesh();

        WATCH(cutMesh->eulerCharacteristic());
        WATCH(cutMesh->isManifold());
        WATCH(cutMesh->nBoundaryLoops());


        int explicitEulerCharacteristic =
            static_cast<int>(static_cast<long long int>(cutMesh->nVertices()) -
                             static_cast<long long int>(cutMesh->nEdges()) +
                             static_cast<long long int>(cutMesh->nFaces()));

        WATCH(explicitEulerCharacteristic);
    }

    ImGui::Separator();
    ImGui::InputText("###CutMeshSaveName", meshSaveName,
                     IM_ARRAYSIZE(meshSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save Cuts")) {
        if (!cutMesh) {
            cutAlongCuts();
        }
        writeSurfaceMesh(*cutMesh, *cutGeometry,
                         std::string(meshSaveName) + ".obj");
    }
}

int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("Conformal cone flattener");
    args::Positional<std::string> meshFilename(
        parser, "mesh", "Mesh to be processed (required).");
    args::ValueFlag<std::string> conesFilename(parser, "string",
                                               "cones filename", {"cones"});
    args::ValueFlag<std::string> curvaturesFilename(
        parser, "string", "curvatures filename", {"curvatures"});
    args::ValueFlag<std::string> scaleFactorsFilename(
        parser, "string", "scaleFactors filename", {"scaleFactors"});
    args::ValueFlag<std::string> ffieldFilename(parser, "string",
                                                "ffield filename", {"ffield"});
    args::ValueFlag<double> greedyConeMaxU(
        parser, "double", "max allowed scale factor in greedy cone placement ",
        {"greedyConeMaxU"});

    args::ValueFlag<std::string> outputMeshFilename(
        parser, "string",
        "obj file to save output mesh to, along with homogeneous texture "
        "coordinates",
        {"outputMeshFilename"});
    args::ValueFlag<std::string> outputLogFilename(
        parser, "string", "csv file to save logs to", {"outputLogFilename"});
    args::Flag useExactCones(parser, "exactCones",
                             "Use exact cones from ffield (no lumping)",
                             {"exactCones"});
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

    // print booleans in words
    std::cout << std::boolalpha;

    if (version) {
        std::cout << "generate cuts version 0.1" << std::endl;
        return 0;
    }

    std::string filename;
    if (meshFilename) {
        filename = args::get(meshFilename);
    } else {
        filename = "../../meshes/thingi10k/nicks_bad_meshes/719791.ply";
        // TODO: better defaults
        // std::cout << "Please provide a mesh file as input." << std::endl;
        // std::cout << parser;
        // return 1;
    }

    // TODO: better defaults
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


    // Load mesh
    std::tie(mesh, geometry) = loadMesh(filename);
    verbose_assert(mesh->nConnectedComponents() == 1, "mesh must be connected");

    std::string nicename = polyscope::guessNiceNameFromPath(filename);

    std::string saveNameGuess = nicename + "_cut";
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
    // Read cones from prescribed curvatures/scale factors if present
    if (conesFilename) {
        loadedCones = true;
        std::cout << "TODO: cone reading not implemented yet" << std::endl;
        // std::pair<size_t, double> prescribedCurvatureInput =
        //     readPrescribedConeAngles(args::get(curvaturesFilename));
        // for (const std::pair<size_t, double>& c :
        // prescribedCurvatureInput)
        //     cones.emplace_back(mesh->vertex(c.first));
    }
    if (curvaturesFilename) {
        loadedCones = true;
        std::vector<std::pair<size_t, double>> prescribedCurvatureInput =
            readPrescribedConeAngles(args::get(curvaturesFilename));
        for (const std::pair<size_t, double>& c : prescribedCurvatureInput)
            cones.emplace_back(mesh->vertex(c.first));
    }
    if (scaleFactorsFilename) {
        loadedCones = true;
        std::vector<std::pair<size_t, double>> prescribedScaleFactorInput =
            readPrescribedScaleFactors(args::get(scaleFactorsFilename));
        for (const std::pair<size_t, double>& c : prescribedScaleFactorInput)
            cones.emplace_back(mesh->vertex(c.first));
    }
    if (ffieldFilename) {
        loadedCones = true;
        std::vector<std::pair<Vertex, double>> ffieldCones =
            readFFieldCones(args::get(ffieldFilename), *mesh, *geometry);

        auto processedCones = args::get(useExactCones)
                                  ? ffieldCones
                                  : lumpCones(*mesh, ffieldCones);

        for (const std::pair<Vertex, double>& c : processedCones)
            cones.emplace_back(c.first);
    }

    if (greedyConeMaxU) {
        uTol = args::get(greedyConeMaxU);
    }

    // TODO: better defaults
    if (args::get(viz) || true) {
        // Initialize polyscope
        polyscope::init();

        // Set the callback function
        polyscope::state::userCallback = myCallback;

        // Register the mesh with polyscope
        psMesh = polyscope::registerSurfaceMesh(
            "input_mesh", geometry->inputVertexPositions,
            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        if (loadedCones) psDrawCones();

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

            bool viz = false;
            generateGoodCuts();
            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        } else {
            std::clock_t start = std::clock();

            bool viz = false;
            generateRandomCones();
            generateGoodCuts();
            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        }

        if (logStats) {
            logger.log("duration", duration);
        }

        if (outputMeshFilename) {
            // TODO: write code to save mesh
        }

        if (outputLogFilename) {
            logger.writeLog(args::get(outputLogFilename));
        }
    }

    return EXIT_SUCCESS;
}
