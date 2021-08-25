#include "SurfaceProjectiveParameterizationQuantity.h"

#include "polyscope/file_helpers.h"
#include "polyscope/polyscope.h"
#include "polyscope/render/engine.h"

#include "imgui.h"

using std::cout;
using std::endl;

namespace polyscope {

// == Custom shader rule for projective texture interpolation
const render::ShaderReplacementRule MESH_PROPAGATE_PROJECTIVE_VALUE2(
    /* rule name */ "MESH_PROPAGATE_PROJECTIVE_VALUE2",
    {
        /* replacement sources */
        {"VERT_DECLARATIONS", R"(
          in vec3 a_value3;
          out vec3 a_value3ToFrag;
        )"},
        {"VERT_ASSIGNMENTS", R"(
          a_value3ToFrag = a_value3;
        )"},
        {"FRAG_DECLARATIONS", R"(
          in vec3 a_value3ToFrag;
        )"},
        {"GENERATE_SHADE_VALUE", R"(
          vec2 shadeValue2 = a_value3ToFrag.xy / a_value3ToFrag.z;
        )"},
    },
    /* uniforms */ {},
    /* attributes */
    {
        {"a_value3", render::DataType::Vector3Float},
    },
    /* textures */ {});

// ==============================================================
// ================  Base Parameterization  =====================
// ==============================================================

SurfaceProjectiveParameterizationQuantity::
    SurfaceProjectiveParameterizationQuantity(std::string name,
                                              ParamCoordsType type_,
                                              ParamVizStyle style_,
                                              SurfaceMesh& mesh_)
    : SurfaceParameterizationQuantity(name, type_, style_, mesh_)

{
    // Register a custom shader rule for projective interpolation
    render::engine->registerShaderRule("MESH_PROPAGATE_PROJECTIVE_VALUE2",
                                       MESH_PROPAGATE_PROJECTIVE_VALUE2);
}


void SurfaceProjectiveParameterizationQuantity::draw() {
    if (!isEnabled()) return;

    if (program == nullptr) {
        createProgram();
    }

    // Set uniforms
    parent.setTransformUniforms(*program);
    setProgramUniforms(*program);
    parent.setStructureUniforms(*program);

    program->draw();
}

void SurfaceProjectiveParameterizationQuantity::createProgram() {
    // Create the program to draw this quantity

    switch (getStyle()) {
    case ParamVizStyle::CHECKER:
        program = render::engine->requestShader(
            "MESH",
            parent.addStructureRules(
                {"MESH_PROPAGATE_PROJECTIVE_VALUE2", "SHADE_CHECKER_VALUE2"}));
        break;
    case ParamVizStyle::GRID:
        program = render::engine->requestShader(
            "MESH",
            parent.addStructureRules(
                {"MESH_PROPAGATE_PROJECTIVE_VALUE2", "SHADE_GRID_VALUE2"}));
        break;
    case ParamVizStyle::LOCAL_CHECK:
        program = render::engine->requestShader(
            "MESH", parent.addStructureRules(
                        {"MESH_PROPAGATE_PROJECTIVE_VALUE2",
                         "SHADE_COLORMAP_ANGULAR2", "CHECKER_VALUE2COLOR"}));
        program->setTextureFromColormap("t_colormap", cMap.get());
        break;
    case ParamVizStyle::LOCAL_RAD:
        program = render::engine->requestShader(
            "MESH",
            parent.addStructureRules(
                {"MESH_PROPAGATE_PROJECTIVE_VALUE2", "SHADE_COLORMAP_ANGULAR2",
                 "SHADEVALUE_MAG_VALUE2", "ISOLINE_STRIPE_VALUECOLOR"}));
        program->setTextureFromColormap("t_colormap", cMap.get());
        break;
    }

    // Fill color buffers
    fillColorBuffers(*program);
    parent.fillGeometryBuffers(*program);

    render::engine->setMaterial(*program, parent.getMaterial());
}

// ==============================================================
// ============  Corner Projective Parameterization  ============
// ==============================================================


SurfaceCornerProjectiveParameterizationQuantity::
    SurfaceCornerProjectiveParameterizationQuantity(
        std::string name, std::vector<glm::vec3> coords_, ParamCoordsType type_,
        ParamVizStyle style_, SurfaceMesh& mesh_)
    : SurfaceProjectiveParameterizationQuantity(name, type_, style_, mesh_),
      coords(std::move(coords_)) {}

std::string SurfaceCornerProjectiveParameterizationQuantity::niceName() {
    return name + " (corner projective parameterization)";
}


void SurfaceCornerProjectiveParameterizationQuantity::fillColorBuffers(
    render::ShaderProgram& p) {
    std::vector<glm::vec3> coordVal;
    coordVal.reserve(3 * parent.nFacesTriangulation());

    size_t cornerCount = 0;
    for (size_t iF = 0; iF < parent.nFaces(); iF++) {
        auto& face = parent.faces[iF];
        size_t D   = face.size();

        // implicitly triangulate from root
        size_t cRoot = cornerCount;
        for (size_t j = 1; (j + 1) < D; j++) {
            size_t cB = cornerCount + j;
            size_t cC = cornerCount + ((j + 1) % D);

            coordVal.push_back(coords[cRoot]);
            coordVal.push_back(coords[cB]);
            coordVal.push_back(coords[cC]);
        }

        cornerCount += D;
    }

    // Store data in buffers
    p.setAttribute("a_value3", coordVal);
}

void SurfaceCornerProjectiveParameterizationQuantity::buildHalfedgeInfoGUI(
    size_t heInd) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();
    ImGui::Text("<%g,%g,%g>", coords[heInd].x, coords[heInd].y,
                coords[heInd].z);
    ImGui::NextColumn();
}

} // namespace polyscope
