#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include "Cutter.h"
#include "Triangulation.h"
#include "Utils.h"

#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"


namespace CEPS {

// Lay out a triangulation in the plane. The non-flat vertices must be provided
// in the cones list. The parameterization is scaled to have an area of
// targetSurfaceArea
CornerData<Vector2> layOutTriangulation(Triangulation& tri,
                                        const std::vector<Vertex>& cones,
                                        double targetSurfaceArea = 1);

// Lay out a triangulation. While cutting the triangulation open, this tries to
// use cuts which have low cost according to edgeCutCosts, although the cuts are
// not guaranteed to be optimal.
CornerData<Vector2> layOutTriangulation(Triangulation& tri,
                                        const std::vector<Vertex>& cones,
                                        const EdgeData<double>& edgeCutCosts,
                                        double targetSurfaceArea = 1);
} // namespace CEPS
