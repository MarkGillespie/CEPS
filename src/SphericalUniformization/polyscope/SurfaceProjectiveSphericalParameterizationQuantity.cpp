#include "SurfaceProjectiveSphericalParameterizationQuantity.h"

#include "polyscope/file_helpers.h"
#include "polyscope/polyscope.h"
#include "polyscope/render/engine.h"

#include <glm/gtx/euler_angles.hpp>

#include "imgui.h"

using std::cout;
using std::endl;

namespace polyscope {

// == Custom shader rule for projective texture interpolation
const render::ShaderReplacementRule MESH_PROPAGATE_PROJECTIVE_VALUE4(
    /* rule name */ "MESH_PROPAGATE_PROJECTIVE_VALUE2",
    {
        /* replacement sources */
        {"VERT_DECLARATIONS", R"(
          in vec4 a_value4;
          out vec4 a_value4ToFrag;
        )"},
        {"VERT_ASSIGNMENTS", R"(
          a_value4ToFrag = a_value4;
        )"},
        {"FRAG_DECLARATIONS", R"(
          in vec4 a_value4ToFrag;
          uniform mat4 u_rot;
          const float PI = 3.1415926535897932384626433832795;
        )"},
        {"GENERATE_SHADE_VALUE", R"(
          vec3 p = vec3(u_rot * vec4(a_value4ToFrag.xyz, 0));
          vec3 projectivePosition = p / a_value4ToFrag.w;
          float theta = atan(projectivePosition.y, projectivePosition.x);
          float phi = acos(projectivePosition.z);
          vec2 tCoord = vec2((theta+PI) / (2*PI), phi / PI);
        )"},
    },
    /* uniforms */ {{"u_rot", render::DataType::Matrix44Float}},
    /* attributes */
    {
        {"a_value4", render::DataType::Vector4Float},
    },
    /* textures */ {});

const render::ShaderReplacementRule SHADE_TEXTURE(
    /* rule name */ "SHADE_TEXTURE",
    {/* replacement sources */
     {"FRAG_DECLARATIONS", R"(
          uniform sampler2D t_ceps_image;
          vec3 undoGammaCorrect(vec3 colorLinear);
        )"},
     {"GENERATE_SHADE_COLOR", R"(
        vec3 albedoColor = undoGammaCorrect(texture(t_ceps_image, tCoord).rgb);
      )"}},
    /* uniforms */
    {},
    /* attributes */ {},
    /* textures */ {{"t_ceps_image", 2}});

// ==============================================================
// =======  Base Projective Spherical Parameterization ==========
// ==============================================================

SurfaceProjectiveSphericalParameterizationQuantity::
    SurfaceProjectiveSphericalParameterizationQuantity(std::string name,
                                                       SurfaceMesh& mesh_)
    : SurfaceMeshQuantity(name, mesh_, true),
      eulerAngles(uniquePrefix() + "#eulerAngles", glm::vec3(0, 0, 0))

{

    // Register a custom shader rule for projective interpolation
    render::engine->registerShaderRule("MESH_PROPAGATE_PROJECTIVE_VALUE4",
                                       MESH_PROPAGATE_PROJECTIVE_VALUE4);
    render::engine->registerShaderRule("SHADE_TEXTURE", SHADE_TEXTURE);


    // Load texture
    int width, height, numComponents;

    unsigned char* imgData =
        stbi_load_from_memory(reinterpret_cast<const unsigned char*>(
                                  &gl::bindata_map_equirectangular[0]),
                              gl::bindata_map_equirectangular.size(), &width,
                              &height, &numComponents, STBI_rgb);
    // unsigned char* imgData =
    //     stbi_load("/Users/mgillesp/Downloads/earth.png", &width, &height,
    //               &numComponents, STBI_rgb);

    if (!imgData) throw std::logic_error("failed to load earth texture image");
    textureBuffer = render::engine->generateTextureBuffer(
        TextureFormat::RGB8, width, height, imgData);
    stbi_image_free(imgData);
}

void SurfaceProjectiveSphericalParameterizationQuantity::draw() {
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

void SurfaceProjectiveSphericalParameterizationQuantity::createProgram() {
    // Create the program to draw this quantity
    program = render::engine->requestShader(
        "MESH", parent.addStructureRules(
                    {"MESH_PROPAGATE_PROJECTIVE_VALUE4", "SHADE_TEXTURE"}));
    // "SHADE_CHECKER_VALUE2"}));
    // Fill color buffers
    fillColorBuffers(*program);
    parent.fillGeometryBuffers(*program);

    render::engine->setMaterial(*program, parent.getMaterial());

    program->setTextureFromBuffer("t_ceps_image", textureBuffer.get());
}


// Update range uniforms
void SurfaceProjectiveSphericalParameterizationQuantity::setProgramUniforms(
    render::ShaderProgram& program) {

    glm::mat4 rotation = glm::eulerAngleYXZ(
        eulerAngles.get().x, eulerAngles.get().y, eulerAngles.get().z);
    program.setUniform("u_rot", &rotation[0][0]);
}

void SurfaceProjectiveSphericalParameterizationQuantity::buildCustomUI() {
    bool updated = false;
    if (ImGui::SliderFloat("alpha", &eulerAngles.get().x, 0.f, 2 * M_PI,
                           "%.3f"))
        updated = true;
    if (ImGui::SliderFloat("beta", &eulerAngles.get().y, 0.f, 2 * M_PI, "%.3f"))
        updated = true;
    if (ImGui::SliderFloat("gamma", &eulerAngles.get().z, 0.f, 2 * M_PI,
                           "%.3f"))
        updated = true;

    if (updated) setEulerAngles(getEulerAngles());
}


SurfaceProjectiveSphericalParameterizationQuantity*
SurfaceProjectiveSphericalParameterizationQuantity::setEulerAngles(
    const glm::vec3& newEulerAngles) {
    eulerAngles = newEulerAngles;
    program.reset();
    requestRedraw();
    return this;
}

glm::vec3 SurfaceProjectiveSphericalParameterizationQuantity::getEulerAngles() {
    return eulerAngles.get();
}
void SurfaceProjectiveSphericalParameterizationQuantity::refresh() {
    program.reset();
    Quantity::refresh();
}

// ==============================================================
// ===============  Vertex Parameterization  ====================
// ==============================================================


SurfaceVertexProjectiveSphericalParameterizationQuantity::
    SurfaceVertexProjectiveSphericalParameterizationQuantity(
        std::string name, std::vector<glm::vec4> coords_, SurfaceMesh& mesh_)
    : SurfaceProjectiveSphericalParameterizationQuantity(name, mesh_),
      coords(std::move(coords_)) {}

std::string
SurfaceVertexProjectiveSphericalParameterizationQuantity::niceName() {
    return name + " (vertex proj-sph parameterization)";
}

void SurfaceVertexProjectiveSphericalParameterizationQuantity::fillColorBuffers(
    render::ShaderProgram& p) {
    std::vector<glm::vec4> coordVal;
    coordVal.reserve(3 * parent.nFacesTriangulation());

    for (size_t iF = 0; iF < parent.nFaces(); iF++) {
        auto& face = parent.faces[iF];
        size_t D   = face.size();

        // implicitly triangulate from root
        size_t vRoot = face[0];
        for (size_t j = 1; (j + 1) < D; j++) {
            size_t vB = face[j];
            size_t vC = face[(j + 1) % D];

            coordVal.push_back(coords[vRoot]);
            coordVal.push_back(coords[vB]);
            coordVal.push_back(coords[vC]);
        }
    }

    // Store data in buffers
    p.setAttribute("a_value4", coordVal);
}

void SurfaceVertexProjectiveSphericalParameterizationQuantity::
    buildVertexInfoGUI(size_t vInd) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();
    ImGui::Text("<%g,%g,%g,%g>", coords[vInd].x, coords[vInd].y, coords[vInd].z,
                coords[vInd].w);
    ImGui::NextColumn();
}


} // namespace polyscope
