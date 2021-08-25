#pragma once

#include "polyscope/affine_remapper.h"
#include "polyscope/histogram.h"
#include "polyscope/render/color_maps.h"
#include "polyscope/render/engine.h"
#include "polyscope/surface_mesh.h"

#include "polyscope/bindata.h"

#include "stb_image.h"

namespace polyscope {


// ==============================================================
// =======  Base Projective Spherical Parameterization ==========
// ==============================================================

class SurfaceProjectiveSphericalParameterizationQuantity
    : public SurfaceMeshQuantity {

  public:
    SurfaceProjectiveSphericalParameterizationQuantity(std::string name,
                                                       SurfaceMesh& mesh_);

    void draw() override;

    virtual void buildCustomUI() override;

    virtual void refresh() override;


    // === Members

    // === Viz stuff
    // to keep things simple, has settings for all of the viz styles, even
    // though not all are used at all times


    // === Getters and setters for visualization options
    SurfaceProjectiveSphericalParameterizationQuantity*
    setEulerAngles(const glm::vec3& newEulerAngles);
    glm::vec3 getEulerAngles();

  protected:
    // === Visualiztion options
    PersistentValue<glm::vec3> eulerAngles;
    std::shared_ptr<render::ShaderProgram> program;
    std::shared_ptr<render::TextureBuffer> textureBuffer;

    // Helpers
    void createProgram();
    void setProgramUniforms(render::ShaderProgram& program);
    virtual void fillColorBuffers(render::ShaderProgram& p) = 0;
};


// ==============================================================
// ===============  Vertex Parameterization  ====================
// ==============================================================

class SurfaceVertexProjectiveSphericalParameterizationQuantity
    : public SurfaceProjectiveSphericalParameterizationQuantity {

  public:
    SurfaceVertexProjectiveSphericalParameterizationQuantity(
        std::string name, std::vector<glm::vec4> values_, SurfaceMesh& mesh_);

    virtual void buildVertexInfoGUI(size_t vInd) override;
    virtual std::string niceName() override;

    // === Members
    std::vector<glm::vec4> coords; // on vertices

  protected:
    virtual void fillColorBuffers(render::ShaderProgram& p) override;
};

template <typename T>
SurfaceVertexProjectiveSphericalParameterizationQuantity*
addProjectiveSphericalParameterizationQuantity(SurfaceMesh& mesh,
                                               std::string name, const T& data);

} // namespace polyscope

#include "SurfaceProjectiveSphericalParameterizationQuantity.ipp"
