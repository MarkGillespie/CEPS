#pragma once

#include "polyscope/affine_remapper.h"
#include "polyscope/histogram.h"
#include "polyscope/render/color_maps.h"
#include "polyscope/render/engine.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_parameterization_enums.h"


namespace polyscope {


// ==============================================================
// ================  Base Parameterization  =====================
// ==============================================================


class SurfaceProjectiveParameterizationQuantity
    : public SurfaceParameterizationQuantity {

  public:
    SurfaceProjectiveParameterizationQuantity(std::string name,
                                              ParamCoordsType type_,
                                              ParamVizStyle style,
                                              SurfaceMesh& mesh_);

    void draw() override;

  protected:
    // Helpers
    void createProgram();
};


// ==============================================================
// ===============  Corner Parameterization  ====================
// ==============================================================

class SurfaceCornerProjectiveParameterizationQuantity
    : public SurfaceProjectiveParameterizationQuantity {

  public:
    SurfaceCornerProjectiveParameterizationQuantity(
        std::string name, std::vector<glm::vec3> values_, ParamCoordsType type_,
        ParamVizStyle style, SurfaceMesh& mesh_);

    virtual void buildHalfedgeInfoGUI(size_t heInd) override;
    virtual std::string niceName() override;

    // === Members
    std::vector<glm::vec3> coords; // on corners

  protected:
    virtual void fillColorBuffers(render::ShaderProgram& p) override;
};

template <typename T>
SurfaceCornerProjectiveParameterizationQuantity*
addProjectiveParameterizationQuantity(
    SurfaceMesh& mesh, std::string name, const T& data,
    ParamCoordsType type = ParamCoordsType::UNIT);

} // namespace polyscope

#include "SurfaceProjectiveParameterizationQuantity.ipp"
